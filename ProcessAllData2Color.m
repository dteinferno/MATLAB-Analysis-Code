%% Clear out the workspace
clear;
clc;

%% Specify the directory
allPathname = uigetdir('C:\Users\turnerevansd\Documents\','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

%% Step through the files, and get an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif')
        tifName
        [RROIs,GROIs] = Get2ColorROIs(allPathname,tifName,-70);
        save(strcat(tifName(1:end-4),'_ROIs','.mat'), 'RROIs', 'GROIs');
        clear RROIs;
        clear GROIs;
    end
end

%% Use the ROIs to get the DF/F for each ROI for the data
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif')
        % Extract the position and DF/F values
        for moveID = 3:length(fileNames)
            moveName =fileNames(moveID).name;
            tifNameParts = strsplit(tifName,'_');
            moveNameParts = strsplit(moveName,'_');
            % Find the associated position data 
            if (strcmpi(moveName(end-3:end),'.txt') &...
                    strcmp(tifName(end-5:end-4),moveName(end-5:end-4)) &...
                    strcmp(tifName(1:4),moveName(1:4)) &...
                    strcmp(moveNameParts{end-1},tifNameParts{end-1}))
                moveName
                [RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] = ImDatLoadBigtiff2Color(tifName,allPathname,0);
                positionDat = VRDatLoad(moveName,allPathname,0);
                clear RstackMean GstackMean;
                load(strcat(tifName(1:end-4),'_ROIs','.mat'));
                num_ROIs=size(RROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(RstackMaxInt);

                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
                GstackXYfiltMax = double(zeros(size(GstackMaxInt)));


                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(RstackMaxInt,3)
                    if mod(i,100)==0
                        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
                    end
                    RstackXYfiltMax(:,:,i) = imfilter(RstackMaxInt(:,:,i),Gxy,'replicate','conv');
                    GstackXYfiltMax(:,:,i) = imfilter(GstackMaxInt(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                RROIaveREFMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                GROIaveREFMax = zeros(num_ROIs,size(GstackXYfiltMax,3));
                h = waitbar(0.0,'Calculating ROIs...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Calculating ROIs...');
                for i = 1:size(RstackXYfiltMax,3)
                    if mod(i,10)==0
                        waitbar(i/size(RstackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltMax,3))]);
                    end
                    for incROI = 1:num_ROIs
                        AR = squeeze(RstackXYfiltMax(:,:,i));
                        AG = squeeze(GstackXYfiltMax(:,:,i));
                        RROIMax = AR(logical(squeeze(RROIs(incROI,:,:))));
                        GROIMax = AG(logical(squeeze(GROIs(incROI,:,:))));
                        RROIaveREFMax(incROI,i) = mean2(RROIMax);
                        GROIaveREFMax(incROI,i) = mean2(GROIMax);
                    end
                end
                delete(h);
                
                RROIaveMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                RROILowMax = RROIaveMax;
                GROIaveMax = zeros(num_ROIs,size(GstackXYfiltMax,3));
                GROILowMax = GROIaveMax;
                for incROI = 1:num_ROIs
                    RROILowMax = sort(squeeze(RROIaveREFMax(incROI,:)));
                    GROILowMax = sort(squeeze(GROIaveREFMax(incROI,:)));
                    RROIaveMax(incROI,:) = RROIaveREFMax(incROI,:)./mean(RROILowMax(1:floor(end/10)));
                    GROIaveMax(incROI,:) = GROIaveREFMax(incROI,:)./mean(GROILowMax(1:floor(end/10)));
                end
      
                save(strcat(tifName(1:end-4),'.mat'), 'RROIaveMax','GROIaveMax','positionDat');
                
                % Savitzky-Golay Filter the RFs
                sgolayOrder = 3;
                sgolayWindow = 11;
                RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
                GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);

                % Plot the profiles
                act = figure('units','normalized','outerposition',[0 0 1 1]);

                % Plot the activity
                subplot(4,1,1);
                imagesc(RROIaveMaxFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Red Maximum Intensity');
                xlabel('Time (sec)');

                subplot(4,1,3);
                imagesc(GROIaveMaxFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Green Maximum Intensity');
                xlabel('Time (sec)');

                colormap(brewermap(64, 'Blues'));

                set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                print(act,strcat(tifName(1:end-4),'_BumpActNONO'),'-dpdf'); 
                
                delete(act);

                clear fullpath RstackMaxInt GstackMaxInt...
                    RstackXYfiltMax GstackXYfiltMax...
                    RROIaveMax GROIaveMax...
                    positionDat;
            end
        end
    end
end    
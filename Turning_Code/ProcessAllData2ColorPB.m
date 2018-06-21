clear;
clc;

% Specify the directory
allPathname = uigetdir('D:/Imaging','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

% Get the reference ROIs
numFlies = input('Number of flies? ');
for fly = 1:numFlies
    PBROIs();
end

% Step through the files, and set an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifName
        [ROIs] = Set2ColorROIsPB(allPathname,tifName);
        save(strcat(tifName(1:end-4),'_ROIs','.mat'), 'ROIs');
        clear ROIs;
    end
end

% Use the ROIs to get the DF/F for each ROI for the data
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifChunks = 1;
        % Find how many tifs are in each stack
%         for tifID2 = 3:length(fileNames)
%             tifName2 = fileNames(tifID2).name;
%             if strcmp(tifName(1:end-5),tifName2(1:end-5))
%                 tifChunks = max(tifChunks,str2num(tifName2(end-4)));
%             end
%         end
        % Extract the position and DF/F values
        for moveID = 3:length(fileNames)
            moveName = fileNames(moveID).name;
            tifNameParts = strsplit(tifName,'_');
            moveNameParts = strsplit(moveName,'_');
            % Find the associated position data 
            if (strcmpi(moveName(end-3:end),'.txt') &...
                    strcmp(tifName(end-5:end-4),moveName(end-5:end-4)) &...
                    strcmp(tifName(1:4),moveName(1:4)) &...
                    strcmp(moveNameParts{end-1},tifNameParts{end-1}))
                tifName
                moveName
                [fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,0);
                load(strcat(tifName(1:end-4),'_ROIs','.mat'));
                num_ROIs=size(ROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(RstackMaxInt);
                
                % Correct for movement
%                 RstackRegMaxInt = imRegSimple(RstackMaxInt, 7);
%                 GstackRegMaxInt = imRegSimple(GstackMaxInt, 7);
%                 RstackRegMean = imRegSimple(RstackMean, 7);
%                 GstackRegMean = imRegSimple(GstackMean, 7);
                RstackRegMaxInt = RstackMaxInt;
                GstackRegMaxInt = GstackMaxInt;
                RstackRegMean = RstackMean;
                GstackRegMean = GstackMean;

                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
                GstackXYfiltMax = double(zeros(size(GstackMaxInt)));

                RstackXYfiltMean = double(zeros(size(RstackMean)));
                GstackXYfiltMean = double(zeros(size(GstackMean)));

                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(RstackMaxInt,3)
                    if mod(i,100)==0
                        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
                    end
                    RstackXYfiltMax(:,:,i) = imfilter(RstackRegMaxInt(:,:,i),Gxy,'replicate','conv');
                    GstackXYfiltMax(:,:,i) = imfilter(GstackRegMaxInt(:,:,i),Gxy,'replicate','conv');
                    RstackXYfiltMean(:,:,i) = imfilter(RstackRegMean(:,:,i),Gxy,'replicate','conv');
                    GstackXYfiltMean(:,:,i) = imfilter(GstackRegMean(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                RROIaveREFMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                GROIaveREFMax = zeros(num_ROIs,size(GstackXYfiltMax,3));
                RROIaveREFMean = zeros(num_ROIs,size(RstackXYfiltMean,3));
                GROIaveREFMean = zeros(num_ROIs,size(GstackXYfiltMean,3));
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
                        BR = squeeze(RstackXYfiltMean(:,:,i));
                        BG = squeeze(GstackXYfiltMean(:,:,i));
                        RROIMax = AR(logical(squeeze(ROIs(incROI,:,:))));
                        GROIMax = AG(logical(squeeze(ROIs(incROI,:,:))));
                        RROIMean = BR(logical(squeeze(ROIs(incROI,:,:))));
                        GROIMean = BG(logical(squeeze(ROIs(incROI,:,:))));
                        RROIaveREFMax(incROI,i) = mean2(RROIMax);
                        GROIaveREFMax(incROI,i) = mean2(GROIMax);
                        RROIaveREFMean(incROI,i) = mean2(RROIMean);
                        GROIaveREFMean(incROI,i) = mean2(GROIMean);
                    end
                end
                delete(h);
                
                RROIaveMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                RROILowMax = RROIaveMax;
                GROIaveMax = zeros(num_ROIs,size(GstackXYfiltMax,3));
                GROILowMax = GROIaveMax;
                RROIaveMean = zeros(num_ROIs,size(RstackXYfiltMean,3));
                RROILowMean = RROIaveMean;
                GROIaveMean = zeros(num_ROIs,size(GstackXYfiltMean,3));
                GROILowMax = GROIaveMean;
                for incROI = 1:num_ROIs
                    RROILowMax = sort(squeeze(RROIaveREFMax(incROI,:)));
                    GROILowMax = sort(squeeze(GROIaveREFMax(incROI,:)));
                    RROILowMean = sort(squeeze(RROIaveREFMean(incROI,:)));
                    GROILowMean = sort(squeeze(GROIaveREFMean(incROI,:)));
                    RROIaveMax(incROI,:) = RROIaveREFMax(incROI,:)./mean(RROILowMax(1:floor(end/10)));
                    GROIaveMax(incROI,:) = GROIaveREFMax(incROI,:)./mean(GROILowMax(1:floor(end/10)));
                    RROIaveMean(incROI,:) = RROIaveREFMean(incROI,:)./mean(RROILowMean(1:floor(end/10)));
                    GROIaveMean(incROI,:) = GROIaveREFMean(incROI,:)./mean(GROILowMean(1:floor(end/10)));
                end
      
                save(strcat(tifName(1:end-4),'.mat'), 'fullpath','RROIaveMax','GROIaveMax','RROIaveMean','GROIaveMean','positionDat');
                
                % Savitzky-Golay Filter the RFs
                sgolayOrder = 3;
                sgolayWindow = 7;
                RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
                GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);
                RROIaveMeanFilt = sgolayfilt(RROIaveMean,sgolayOrder,sgolayWindow,[],2);
                GROIaveMeanFilt = sgolayfilt(GROIaveMean,sgolayOrder,sgolayWindow,[],2);

                % Plot the profiles
                act = figure('units','normalized','outerposition',[0 0 1 1]);

                % Plot the activity
                subplot(4,1,1);
                imagesc(RROIaveMaxFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Red Maximum Intensity');
                xlabel('Time (sec)');

                subplot(4,1,2);
                imagesc(RROIaveMeanFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Red Mean Intensity');
                xlabel('Time (sec)');

                subplot(4,1,3);
                imagesc(GROIaveMaxFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Green Maximum Intensity');
                xlabel('Time (sec)');

                subplot(4,1,4);
                imagesc(GROIaveMeanFilt);
                hold on;
                set(gca,'FontSize',16);
                title('Green Mean Intensity');
                xlabel('Time (sec)');

                colormap(brewermap(64, 'Blues'));

                set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                print(act,strcat(tifName(1:end-4),'_BumpAct'),'-dpdf'); 
                
                delete(act);

                clear fullpath RstackMean GstackMean RstackMaxInt GstackMaxInt positionDat...
                    RROIaveMax GROIaveMax RROIaveMean GROIaveMean;
            end
        end
    end
end    
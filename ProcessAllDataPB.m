clear;
clc;

% Specify the directory
allPathname = uigetdir('C:\Users\turnerevansd\Documents','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

% Get the reference ROIs
numFlies = input('Number of flies? ');
for fly = 1:numFlies
    PBROIs(1);
end

%% Step through the files, and set an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifName
        [ROIs] = SetROIsPB(allPathname,tifName);
        save(strcat(tifName(1:end-4),'_ROIs','.mat'), 'ROIs');
        clear ROIs;
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
                [stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathname,0);
                clear stackMean;
                stackReg = imRegSimple(stackMaxInt, 10); % Image correct the stacks
                
                positionDat = VRDatLoad(moveName,allPathname,0);
                load(strcat(tifName(1:end-4),'_ROIs','.mat'));
                num_ROIs=size(ROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(stackMaxInt);

                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                stackXYfiltMax = double(zeros(size(stackReg)));

                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(stackReg,3)
                    if mod(i,100)==0
                        waitbar(i/length(stackReg),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stackReg))]);
                    end
                    stackXYfiltMax(:,:,i) = imfilter(stackReg(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                ROIaveREFMax = zeros(num_ROIs,size(stackXYfiltMax,3));

                h = waitbar(0.0,'Calculating ROIs...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Calculating ROIs...');
                for i = 1:size(stackXYfiltMax,3)
                    if mod(i,10)==0
                        waitbar(i/size(stackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(stackXYfiltMax,3))]);
                    end
                    for incROI = 1:num_ROIs
                        A = squeeze(stackXYfiltMax(:,:,i));
                        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
                        ROIaveREFMax(incROI,i) = mean2(ROIMax);
                    end
                end
                delete(h);
                
                ROIaveMax = zeros(num_ROIs,size(stackXYfiltMax,3));
                ROILowMax = ROIaveMax;
                
                for incROI = 1:num_ROIs
                    ROILowMax = sort(squeeze(ROIaveREFMax(incROI,:)));
                    ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./mean(ROILowMax(1:floor(end/10)));
                end
                
                save(strcat(tifName(1:end-4),'.mat'), 'ROIaveMax','positionDat');
                
                % Savitzky-Golay Filter the RFs
                sgolayOrder = 3;
                sgolayWindow = 11;
                ROIaveMaxFilt = sgolayfilt(ROIaveMax,sgolayOrder,sgolayWindow,[],2);
                
                % Plot the profiles
                act = figure('units','normalized','outerposition',[0 0 1 1]);

                % Plot the activity
                subplot(4,1,1);
                imagesc(ROIaveMaxFilt);
                hold on;
                set(gca,'FontSize',16);
                title('maximum intensity');
                xlabel('time (sec)');

                colormap(brewermap(64, 'Blues'));

                set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                print(act,strcat(tifName(1:end-4),'_BumpAct'),'-dpdf'); 
                
                delete(act);

                clear fullpath stackMaxInt stackReg stackXYfiltMax ROIaveMax positionDat;
            end
        end
    end
end    
% Correct NO data

%% List directories
dirID{1} = 'D:\Imaging\2Color\20160504';

%% remake ROIs for all NO data FOR TRIALS IN THE DARK

for dirNum = 1
    % Specify the directory
    allPathname = dirID{dirNum};
    fileNames = dir(allPathname);
    cd(allPathname);
    allPathname = strcat(allPathname,'\');

    % Step through the files, and get an ROI for each set of tifs
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-4:end),'1.tif') & strcmp(tifName(end-19:end-16),'Dark')
            tifName
            [RROIs,GROIs] = Get2ColorROIs(allPathname,tifName);
            save(strcat(tifName(1:end-10),'_ROIs','.mat'), 'RROIs','-append');
            clear RROIs;
        end
    end

end

%% Rerun ROI calc for all NO data FOR TRIALS IN THE DARK

for dirNum = 1
    % Specify the directory
    allPathname = dirID{dirNum};
    fileNames = dir(allPathname);
    cd(allPathname);
    allPathname = strcat(allPathname,'\');

% Use the ROIs to get the DF/F for each ROI for the data
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-4:end),'1.tif') & strcmp(tifName(end-19:end-16),'Dark')
        tifChunks = 1;
        % Find how many tifs are in each stack
        for tifID2 = 3:length(fileNames)
            tifName2 = fileNames(tifID2).name;
            if strcmp(tifName(1:end-5),tifName2(1:end-5))
                tifChunks = max(tifChunks,str2num(tifName2(end-4)));
            end
        end
        % Extract the position and DF/F values
        for moveID = 3:length(fileNames)
            moveName =fileNames(moveID).name;
            tifNameParts = strsplit(tifName,'_');
            moveNameParts = strsplit(moveName,'_');
            % Find the associated position data 
            if (strcmpi(moveName(end-3:end),'.txt') &...
                    strcmp(tifName(end-11:end-10),moveName(end-5:end-4)) &...
                    strcmp(tifName(1:4),moveName(1:4)) &...
                    strcmp(moveNameParts{end-1},tifNameParts{end-2}))
                [fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,1);
                load(strcat(tifName(1:end-10),'_ROIs','.mat'));
                num_ROIs=size(RROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(RstackMaxInt);
                
                % Correct for movement
                RstackRegMaxInt = imRegSimple(RstackMaxInt, 20);
                RstackRegMean = imRegSimple(RstackMean, 20);

                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
                RstackXYfiltMean = double(zeros(size(RstackMean)));

                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(RstackMaxInt,3)
                    if mod(i,100)==0
                        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
                    end
                    RstackXYfiltMax(:,:,i) = imfilter(RstackRegMaxInt(:,:,i),Gxy,'replicate','conv');
                    RstackXYfiltMean(:,:,i) = imfilter(RstackRegMean(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                RROIaveREFMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                RROIaveREFMean = zeros(num_ROIs,size(RstackXYfiltMean,3));
                h = waitbar(0.0,'Calculating ROIs...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Calculating ROIs...');
                for i = 1:size(RstackXYfiltMax,3)
                    if mod(i,10)==0
                        waitbar(i/size(RstackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltMax,3))]);
                    end
                    for incROI = 1:num_ROIs
                        AR = squeeze(RstackXYfiltMax(:,:,i));
                        BR = squeeze(RstackXYfiltMean(:,:,i));
                        RROIMax = AR(logical(squeeze(RROIs(incROI,:,:))));
                        RROIMean = BR(logical(squeeze(RROIs(incROI,:,:))));
                        RROIaveREFMax(incROI,i) = mean2(RROIMax);
                        RROIaveREFMean(incROI,i) = mean2(RROIMean);
                    end
                end
                delete(h);
                
                RROIaveMax = zeros(num_ROIs,size(RstackXYfiltMax,3));
                RROILowMax = RROIaveMax;
                RROIaveMean = zeros(num_ROIs,size(RstackXYfiltMean,3));
                RROILowMax = RROIaveMean;
                for incROI = 1:num_ROIs
                    RROILowMax = sort(squeeze(RROIaveREFMax(incROI,:)));
                    RROILowMean = sort(squeeze(RROIaveREFMean(incROI,:)));
                    RROIaveMax(incROI,:) = RROIaveREFMax(incROI,:)./mean(RROILowMax(1:floor(end/10)));
                    RROIaveMean(incROI,:) = RROIaveREFMean(incROI,:)./mean(RROILowMean(1:floor(end/10)));
                end
      
                save(strcat(tifName(1:end-10),'.mat'), 'RROIaveMax','RROIaveMean','-append');
                
                clear fullpath RstackMean RstackMaxInt ...
                    RROIaveMax RROIaveMean;
            end
        end
    end
end
end
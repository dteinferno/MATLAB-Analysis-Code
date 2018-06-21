%% Clear out the workspace
clear;
clc;

%% Get the directories and, for each directory, 
% set the number of flies, the region that is being imaged, and whether one or two color imaging was performed.
% Go to the root directory
dirRoot = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the root directory');
cd(dirRoot);
% Specify the directories
numDirs = input('Number of directories? ');
allPathname = cell(numDirs,1);
for dirNow = 1:numDirs
    allPathname{dirNow}.name = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the directory');
    allPathname{dirNow}.numFlies = input('Number of flies? ');
    allPathname{dirNow}.numColors = input('1 or 2 color imaging? ');
    for colorID = 1:allPathname{dirNow}.numColors
        if colorID == 1
            colorName = 'Green';
        else
            colorName = 'Red';
        end
        allPathname{dirNow}.imRegion{colorID} = input(strcat(colorName, ' region that was imaged? (EB,PB,FB,other) '));
        while ~(strcmp(allPathname{dirNow}.imRegion{colorID},'EB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'PB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'FB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'other'))
            allPathname{dirNow}.imRegion{colorID} = input('Region that was imaged? (EB,PB,other)');
        end
    end
end

%% Get the reference ROIs
for dirNow = 1:numDirs
    cd(allPathname{dirNow}.name);
    for fly = 1:allPathname{dirNow}.numFlies
        GetRefROIs(allPathname{dirNow}.imRegion,allPathname{dirNow}.numColors);
    end
end

%% Set the ROIs
for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif') & ~contains(tifName,'Anatomy')
            tifName
            ROIs = SetROIs(allPathnameNow,tifName,allPathname{dirNow}.imRegion,allPathname{dirNow}.numColors);
            
            newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
            ROIfName = strcat(newPath,tifName(1:end-4));
    
            save(strcat(ROIfName,'_ROIs','.mat'), 'ROIs');
            clear ROIs;
        end
    end 
end

%% For each stack in the directory, find the cutoff ROIs 

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif') & ~contains(tifName,'Anatomy')
            % Get the tiff header info
            fullpath = strcat(allPathnameNow,tifName);
            reader = ScanImageTiffReader(fullpath);
            desc = reader.metadata;
            planeLoc = strfind(desc, 'numFramesPerVolume');
            discardLoc = strfind(desc, 'numDiscardFlybackFrames');
            num_planes = str2num(desc(planeLoc+21:planeLoc+22));
            num_discards = str2num(desc(discardLoc+26));

            % Load the tifStack
            vol=ScanImageTiffReader(fullpath).data();
            width = size(vol,1);
            height = size(vol,2);
            num_images = size(vol,3);
            numVolumes = floor(num_images/num_planes);

            % Find the average for each plane
            planeAverage = int16(zeros(num_planes,height,width));
            h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
            set(h,'Position',[50 50 360 72]);
            set(h,'Name','Loading TIFF stack...');
            for incIm = 1:num_images
                if mod(incIm,100)==0
                    waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
                end
                planeAverage(1+mod(incIm-1,num_planes),:,:) = squeeze(planeAverage(1+mod(incIm-1,num_planes),:,:)) + vol(:,:,incIm)';
            end
            delete(h);
            planeAverage = planeAverage./numVolumes;

            % For each plane, specify the region of interest
            cutROIs = zeros(num_planes,height,width);
            for planeNow = 1:num_planes-num_discards
                planePic = figure;
                A = squeeze(planeAverage(planeNow,:,:));
                imagesc(A);
                axis equal;
                cutOut = impoly;
                wait(cutOut);
                position = getPosition(cutOut);
                cutROIs(planeNow,:,:) = roipoly(A,position(:,1),position(:,2));
                delete(planePic);
            end
            newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
            ROIfName = strcat(newPath,tifName(1:end-4));
            save(strcat(ROIfName,'_cutROIs','.mat'), 'cutROIs');
        end
    end
end

%% Use the cut ROIs to restrict the MIP and the region ROIs to get the DF/F

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
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
                        strcmp(moveNameParts{end-1},tifNameParts{end-1}) & ...
                        strcmp(moveNameParts{end-2},tifNameParts{end-2}))
                    moveName
                    
                    % Load the stacks
                    fullpath = strcat(allPathnameNow,tifName);
                    reader = ScanImageTiffReader(fullpath);
                    desc = reader.metadata;
                    planeLoc = strfind(desc, 'numFramesPerVolume');
                    discardLoc = strfind(desc, 'numDiscardFlybackFrames');
                    num_planes = str2num(desc(planeLoc+21:planeLoc+22));
                    num_discards = str2num(desc(discardLoc+26));

                    % Load the tifStack
                    vol=ScanImageTiffReader(fullpath).data();
                    width = size(vol,1);
                    height = size(vol,2);
                    num_images = size(vol,3);
                    numVolumes = floor(num_images/num_planes);
                    
                    % Load the cutoff ROIs
                    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
                    ROIfName = strcat(newPath,tifName(1:end-4));
                    load(strcat(ROIfName,'_cutROIs','.mat'));
                    
                    % Calculate the maximum intensity projection using the
                    % cutoff ROIs
                    stackMaxInt = double(zeros(height,width,numVolumes));
                    miniStack = double(zeros(height,width,(num_planes-num_discards)));
                    h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Loading TIFF stack...');
                    for incIm = 1:num_images
                        if mod(incIm,100)==0
                            waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
                        end
                        if mod(ceil(incIm/2),num_planes) ~= 0
                            miniStack(:,:,mod(ceil(incIm/2),num_planes)) = double(squeeze(vol(:,:,incIm)')).*squeeze(cutROIs(1+mod(incIm-1,num_planes),:,:));
                            if (mod(incIm,num_planes) == num_planes-num_discards)
                                stackMaxInt(:,:,ceil(incIm/(num_planes))) = max(miniStack,[],3);
                            end
                        end
                    end
                    delete(h);
                    
                    if strcmp(allPathname{dirNow}.imRegion{1},'EB') || strcmp(allPathname{dirNow}.imRegion{1},'FB')
                        rotAng = -70;
                    else
                        rotAng = 90;
                    end
                    
                    stackReg = stackMaxInt; % Image correct the stacks
                    % Rotate the stacks
                    stkMean = mean(stackReg,3);
                    stkMean = imrotate(stkMean,rotAng);
                    stackMaxIntRot = zeros([size(stkMean) size(stackReg,3)]);
                    for frame = 1:size(stackReg,3)
                        stackMaxIntRot(:,:,frame) = imrotate(stackReg(:,:,frame),rotAng);
                    end

                    % Load the movement data
                    positionDat = VRDatLoad(moveName,allPathnameNow,0);
                    
                    % Load the ROIs
                    load(strcat(ROIfName,'_ROIs','.mat'));
                    numROIsG=size(ROIs{1},1);
                    numPlanes = length(positionDat.tFrameGrab)/length(stackMaxInt);

                    % Gaussian Filter
                    gaussianSize = [2 2];
                    gaussianSigma = 0.5;
                    Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
                    h = waitbar(0.0,'Gaussian filtering stack...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Gaussian filtering TIFF stack...');

                    stackXYfiltMax = double(zeros(size(stackMaxIntRot)));
                    for i = 1:size(stackMaxIntRot,3)
                        if mod(i,100)==0
                            waitbar(i/length(stackMaxIntRot),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stackMaxIntRot))]);
                        end
                        stackXYfiltMax(:,:,i) = imfilter(stackMaxIntRot(:,:,i),Gxy,'replicate','conv');
                    end
                    delete(h);
                        
                    % Plot the ROI on each figure and calculate the average
                    h = waitbar(0.0,'Calculating ROIs...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Calculating ROIs...');
                    ROIaveREFMax = zeros(numROIsG,size(stackXYfiltMax,3));
                    ROIsNow = squeeze(ROIs{1});

                    for i = 1:size(stackXYfiltMax,3)
                        if mod(i,10)==0
                            waitbar(i/size(stackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(stackXYfiltMax,3))]);
                        end
                        for incROI = 1:numROIsG
                            A = squeeze(stackXYfiltMax(:,:,i));
                            ROIMax = A(logical(squeeze(ROIsNow(incROI,:,:))));
                            ROIaveREFMax(incROI,i) = mean2(ROIMax);
                        end
                    end


                    ROIaveMax = zeros(numROIsG,size(stackXYfiltMax,3));
                    ROILowMax = ROIaveMax;
                    numLow = floor(length(ROILowMax)/10);
                    for incROI = 1:numROIsG
                        ROILowMax = sort(squeeze(ROIaveREFMax(incROI,:)));
                        ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./mean(ROILowMax(1:numLow));
                    end

                    save(strcat(ROIfName,'.mat'), 'ROIaveMax','positionDat');
                    delete(h);

                    % Savitzky-Golay Filter the RFs
                    sgolayOrder = 3;
                    sgolayWindow = 11;
                    ROIaveMaxFilt = sgolayfilt(ROIaveMax,sgolayOrder,sgolayWindow,[],2);

                    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
                    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/numPlanes);
                    maxFG = round(length(positionDat.tFrameGrab)/numPlanes);
                    positionDat.num_planes = numPlanes;
                    positionDat.minFG = minFG;
                    positionDat.maxFG = maxFG;

                    % convert the raw ball position data to position (in case
                    % the visual display is not closed loop with a gain of 1)
                    [posRot, posFor, posLat] = PositionConverter(positionDat);

                    % Match the position data to the framegrab times
                    OffsetRotMatch = MatchData(positionDat.t,positionDat);
                    OffsetRotMatch(:,2) = MatchData(pi/180*positionDat.OffsetRot,positionDat);
                    OffsetForMatch = MatchData(positionDat.OffsetFor,positionDat);
                    OffsetLatMatch = MatchData(positionDat.OffsetLat,positionDat);
                    PosRotMatch = MatchData(posRot,positionDat);
                    PosForMatch = MatchData(posFor,positionDat);
                    PosLatMatch = MatchData(posLat,positionDat);

                    % Plot the profiles
                    act = figure('units','normalized','outerposition',[0 0 1 1]);

                    % Plot the activity
                    subplot(4,1,1);
                    hold on;
                    plot(OffsetRotMatch(:,1),OffsetRotMatch(:,2),'b');
                    plot(OffsetRotMatch(:,1),pi/180*PosRotMatch,'k');
                    legend({'stripe','heading'});
                    ylabel('position (rad)');
                    ylim([-pi pi]);
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                    
                    subplot(4,1,2);
                    imagesc(OffsetRotMatch(:,1),[1:numROIsG],flipud(ROIaveMaxFilt));
                    title('maximum intensity');
                    xlabel('time (sec)');
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                    colormap(brewermap(64, 'Blues'));

                    set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                    print(act,strcat(ROIfName,'_BumpAct'),'-dpdf'); 

                    delete(act);

                    clear fullpath stackMaxInt stackReg stackXYfiltMax ROIaveMax positionDat;
                end
            end
        end
    end    
end
%% Clear out the workspace
clear;
clc;

%% Get the directories and, for each directory, 
% the number of flies, the region that is being imaged, and whether one or two color imaging was performed.
% Go to the root directory
addpath('/Users/loaner/Documents/MATLAB/analysis_code/ScanImageTiffReader-1.1-Darwin/share/matlab');
addpath('/Users/loaner/Documents/MATLAB/analysis_code/ImagingDataScripts')%Kris
addpath('/Users/loaner/Documents/MATLAB/analysis_code/BehaviorScripts')
addpath('/Users/loaner/Documents/MATLAB/analysis_code/DrosteEffect-BrewerMap-04533de');
addpath('/Users/loaner/Documents/MATLAB/analysis_code/DrosteEffect-BrewerMap-04533de');
addpath( '/Users/loaner/Documents/MATLAB/line_specific');
addpath( '/Users/loaner/Documents/MATLAB/analysis_code');
addpath( '/Users/loaner/Documents/MATLAB/analysis_code/distributionPlot/distributionPlot');
%dirRoot = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the root directory');
dirRoot = '/Users/loaner/Documents/MATLAB/analysis_code';
cd(dirRoot);
% Specify the directories
%numDirs = input('Number of directories? ');
numDirs = 1; %Kris
allPathname = cell(numDirs,1);

for dirNow = 1:numDirs
    %allPathname{dirNow}.name = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the directory');
    allPathname{dirNow}.name = '/Users/loaner/Documents/Imaging/G_PEN1-R_EPG_test'; %Kris
    %allPathname{dirNow}.numFlies = input('Number of flies? ');
    allPathname{dirNow}.numFlies = 1; %Kris
    %allPathname{dirNow}.numColors = input('1 or 2 color imaging? ');
    allPathname{dirNow}.numColors = 2; %Kris
    for colorID = 1:allPathname{dirNow}.numColors
        if colorID == 1
            colorName = 'Green';
        else
            colorName = 'Red';
        end
        %allPathname{dirNow}.imRegion{colorID} = input(strcat(colorName, ' region that was imaged? (EB,PB,FB,NO,other) '));
        allPathname{dirNow}.imRegion{colorID} = 'EB'; %Kris
        while ~(strcmp(allPathname{dirNow}.imRegion{colorID},'EB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'PB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'FB') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'NO') ||...
                strcmp(allPathname{dirNow}.imRegion{colorID},'other'))
            allPathname{dirNow}.imRegion{colorID} = input('Region that was imaged? (EB,PB,NO,other)');
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
    allPathnameNow = strcat('/',allPathname{dirNow}.name,'/');
    for tifID = 3:length(fileNames) %for each tif
        tifName = fileNames(tifID).name;
        % Process each trial
        %tifName
        if strcmp(tifName(end-3:end),'.tif') %& ~contains(tifName,'Anatomy') if it's actually a tif
            display(tifName)
            display(allPathnameNow)
            display('new ROI')
            %check ROIs
            ROIs = SetROIs(allPathnameNow,tifName,allPathname{dirNow}.imRegion,allPathname{dirNow}.numColors);
            
            newPath = strrep(allPathnameNow,'Downloads/Data','Documents/RawAnalysis');
            ROIfName = strcat(newPath,tifName(1:end-4));
    
            save(strcat(ROIfName,'_ROIs','.mat'), 'ROIs'); %what are these mats?
            display('saved ROIs')
            clear ROIs;
        end
    end 
end

%% Use the ROIs to get the DF/F for each ROI for the data

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'/');
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        display(tifName)
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif')
            % Extract the position and DF/F values
            for moveID = 3:length(fileNames)

                moveName = fileNames(moveID).name;
                %display(moveName)
                tifNameParts = strsplit(tifName,'_');
                moveNameParts = strsplit(moveName,'_');
                % Find the associated position data 
                if (strcmpi(moveName(end-3:end),'.txt') &...
                        strcmp(tifName(end-5:end-4),moveName(end-5:end-4)) &...
                        strcmp(tifName(1:4),moveName(1:4)) &...
                        strcmp(moveNameParts{end-1},tifNameParts{end-1}) & ...
                        strcmp(moveNameParts{end-2},tifNameParts{end-2}))
                    display('found file')
                    display(moveName)
                    
                    % Load the stacks
                    display('loading stacks')
                    if strcmp(allPathname{dirNow}.imRegion{1},'EB') || strcmp(allPathname{dirNow}.imRegion{1},'FB')
                        rotAng = -70;
                    else
                        rotAng = 90;
                    end
    
                    if allPathname{dirNow}.numColors == 1
                        [stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,0);
                        clear stackMean;
                        stackReg = imRegSimple(stackMaxInt, 10); % Image correct the stacks
                        % Rotate the stacks
                        stkMean = mean(stackReg,3);
                        stkMean = imrotate(stkMean,rotAng);
                        stackMaxIntRot = zeros([size(stkMean) size(stackReg,3)]);
                        for frame = 1:size(stackReg,3)
                            stackMaxIntRot(:,:,frame) = imrotate(stackReg(:,:,frame),rotAng);
                        end
                    elseif allPathname{dirNow}.numColors == 2
                        [RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] =...
                            ImDatLoadBigtiff2Color(tifName,allPathnameNow,0); %256x256 image for each timepoint
                        %this is the data from a max int projection as a
                        %256*256 grid. Calculate normalized RMSD between
                        %these two vs. time/velocity?
                        clear RstackMean GstackMean;
                        [GstackReg RstackReg] = imRegSimple2Color(GstackMaxInt, RstackMaxInt, 10); % Image correct the stacks
                        %these are the image corrected stacks??
                        % Rotate the stacks so that the EB is aligned
                        RMean = mean(RstackReg,3);
                        RMean = imrotate(RMean,rotAng);
                        GMean = mean(GstackReg,3);
                        GMean = imrotate(GMean,rotAng);
                        RstackMaxIntRot = zeros([size(RMean) size(RstackReg,3)]);
                        GstackMaxIntRot = zeros([size(GMean) size(GstackReg,3)]);
                        for frame = 1:size(RstackMaxIntRot,3)
                            RstackMaxIntRot(:,:,frame) = imrotate(RstackReg(:,:,frame),rotAng);
                            GstackMaxIntRot(:,:,frame) = imrotate(GstackReg(:,:,frame),rotAng);
                        end
                        %these are rotated appropriately - finished messing
                        %with them
                    end
                    
                    % Load the movement data
                    positionDat = VRDatLoad(moveName,allPathnameNow,0);
                    
                    display('loaded')
                    
                    % Load the ROIs
                    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
                    ROIfName = strcat(newPath,tifName(1:end-4));
                    load(strcat(ROIfName,'_ROIs','.mat'));
                    numROIsG=size(ROIs{1},1);
                    display('loaded things for Gaussian')
                    if allPathname{dirNow}.numColors == 2
                        numROIsR=size(ROIs{2},1);
                    end
                    display('loaded ROI')
                    if allPathname{dirNow}.numColors == 1
                        numPlanes = length(positionDat.tFrameGrab)/length(stackMaxInt);
                    elseif allPathname{dirNow}.numColors == 2
                        numPlanes = length(positionDat.tFrameGrab)/length(RstackMaxInt);
                    end
                    
                    display('starting Gaussian')

                    % Gaussian Filter
                    gaussianSize = [2 2];
                    gaussianSigma = 0.5;
                    Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
                    h = waitbar(0.0,'Gaussian filtering stack...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Gaussian filtering TIFF stack...');

                    if allPathname{dirNow}.numColors == 1
                        stackXYfiltMax = double(zeros(size(stackMaxIntRot)));
                        for i = 1:size(stackMaxIntRot,3)
                            if mod(i,100)==0
                                waitbar(i/length(stackMaxIntRot),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stackMaxIntRot))]);
                            end
                            stackXYfiltMax(:,:,i) = imfilter(stackMaxIntRot(:,:,i),Gxy,'replicate','conv');
                        end
                    elseif allPathname{dirNow}.numColors == 2
                        RstackXYfiltMax = double(zeros(size(RstackMaxIntRot)));
                        GstackXYfiltMax = double(zeros(size(GstackMaxIntRot)));
                        for i = 1:size(RstackMaxIntRot,3)
                            if mod(i,100)==0
                                waitbar(i/length(RstackMaxIntRot),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxIntRot))]);
                            end
                            RstackXYfiltMax(:,:,i) = imfilter(RstackMaxIntRot(:,:,i),Gxy,'replicate','conv');
                            GstackXYfiltMax(:,:,i) = imfilter(GstackMaxIntRot(:,:,i),Gxy,'replicate','conv');
                        end
                    end
                    delete(h);
                        
                    % Plot the ROI on each figure and calculate the average
                    h = waitbar(0.0,'Calculating ROIs...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Calculating ROIs...');
                    if allPathname{dirNow}.numColors == 1
                        ROIaveREFMax = zeros(numROIsG,size(stackXYfiltMax,3));
                        ROIsNow = squeeze(ROIs{1});
                        
                        for i = 1:size(stackXYfiltMax,3)
                            if mod(i,10)==0
                                waitbar(i/size(stackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(stackXYfiltMax,3))]);
                            end
                            for incROI = 1:numROIsG
                                A = squeeze(stackXYfiltMax(:,:,i));
                                if numROIsG > 1
                                    ROIMax = A(logical(squeeze(ROIsNow(incROI,:,:))));
                                    ROIaveREFMax(incROI,i) = mean2(ROIMax);
                                else
                                    ROIMax = A(logical(ROIsNow));
                                    ROIaveREFMax(incROI,i) = mean2(ROIMax);
                                end
                            end
                        end
                        

                        ROIaveMax = zeros(numROIsG,size(stackXYfiltMax,3));
                        ROILowMax = ROIaveMax;
                        for incROI = 1:numROIsG
                            ROILowMax = sort(squeeze(ROIaveREFMax(incROI,:)));
                            ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./mean(ROILowMax(1:floor(end/10)));
                        end

                        save(strcat(ROIfName,'.mat'), 'ROIaveMax','positionDat');
                    elseif allPathname{dirNow}.numColors == 2
                        RROIaveREFMax = zeros(numROIsR,size(RstackXYfiltMax,3));
                        GROIaveREFMax = zeros(numROIsG,size(GstackXYfiltMax,3));
                        RROIs = squeeze(ROIs{2}); %defines ROIs for red
                        GROIs = squeeze(ROIs{1}); %defines ROIs for green
                        for i = 1:size(RstackXYfiltMax,3)
                            if mod(i,10)==0
                                waitbar(i/size(RstackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltMax,3))]);
                            end
                            for incROI = 1:numROIsG
                                AG = squeeze(GstackXYfiltMax(:,:,i));
                                GROIMax = AG(logical(squeeze(GROIs(incROI,:,:))));
                                GROIaveREFMax(incROI,i) = mean2(GROIMax);
                            end
                            for incROI = 1:numROIsR
                                AR = squeeze(RstackXYfiltMax(:,:,i));
                                RROIMax = AR(logical(squeeze(RROIs(incROI,:,:))));
                                RROIaveREFMax(incROI,i) = mean2(RROIMax);
                            end
                        end

                        RROIaveMax = zeros(numROIsR,size(RstackXYfiltMax,3));
                        RROILowMax = RROIaveMax;
                        GROIaveMax = zeros(numROIsG,size(GstackXYfiltMax,3));
                        GROILowMax = GROIaveMax;
                        for incROI = 1:numROIsG
                            GROILowMax = sort(squeeze(GROIaveREFMax(incROI,:)));
                            GROIaveMax(incROI,:) = GROIaveREFMax(incROI,:)./mean(GROILowMax(1:floor(end/10)));
                        end
                        for incROI = 1:numROIsR
                            RROILowMax = sort(squeeze(RROIaveREFMax(incROI,:)));
                            RROIaveMax(incROI,:) = RROIaveREFMax(incROI,:)./mean(RROILowMax(1:floor(end/10)));
                            %average max intensity over the ROI??
                        end

                        save(strcat(ROIfName,'.mat'), 'RROIaveMax','GROIaveMax','positionDat');
                    end
                    delete(h);

                    % Savitzky-Golay Filter the RFs
                    sgolayOrder = 3;
                    sgolayWindow = 11;
                    if allPathname{dirNow}.numColors == 1
                        ROIaveMaxFilt = sgolayfilt(ROIaveMax,sgolayOrder,sgolayWindow,[],2);
                    elseif allPathname{dirNow}.numColors == 2
                        RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
                        GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);
                    end
                    
                    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000; %corrected times
                    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/numPlanes);
                    maxFG = round(length(positionDat.tFrameGrab)/numPlanes);
                    %divide by number of planes to get timepoints rather
                    %than frames (e.g. 10 frames per cycle - divide by 10)
                    positionDat.num_planes = numPlanes;
                    positionDat.minFG = minFG;%time of start and stop of imaging?
                    positionDat.maxFG = maxFG;

                    % convert the raw ball position data to position (in case
                    % the visual display is not closed loop with a gain of 1)
                    [posRot, posFor, posLat] = PositionConverter(positionDat); %what are these? posRot is fly orientation
                    %postition as a function of time but how is it quantified?
                    
                    display('got positions')

                    % Match the position data to the framegrab times
                    OffsetRotMatch = MatchData(positionDat.t,positionDat); %store times
                    OffsetRotMatch(:,2) = MatchData(pi/180*positionDat.OffsetRot,positionDat); %store stripe position
                    OffsetForMatch = MatchData(positionDat.OffsetFor,positionDat);
                    OffsetLatMatch = MatchData(positionDat.OffsetLat,positionDat);
                    PosRotMatch = MatchData(posRot,positionDat); %adjusts times - takes data at times
                    PosForMatch = MatchData(posFor,positionDat);
                    PosLatMatch = MatchData(posLat,positionDat);

                    % Plot the profiles
                    act = figure('units','normalized','outerposition',[0 0 1 1]);

                    % Plot the activity
                    subplot(4,1,1);
                    hold on;
                    plot(OffsetRotMatch(:,1),OffsetRotMatch(:,2),'b'); %offset vs time. Offset is position of stripe
                    plot(OffsetRotMatch(:,1),pi/180*PosRotMatch,'k'); %this is the actual heading
                    legend({'stripe','heading'});
                    ylabel('position (rad)');
                    ylim([-pi pi]);
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                    
                    
                    if allPathname{dirNow}.numColors == 1
                        subplot(4,1,2);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsG],flipud(ROIaveMaxFilt));
                        %time, range(ROIs), flip ROIs upsidedown.
                        title('maximum intensity');
                        xlabel('time (sec)');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                    elseif allPathname{dirNow}.numColors == 2
                        subplot(4,1,2);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsR],flipud(RROIaveMaxFilt));
                        title('red maximum intensity');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
               
                        subplot(4,1,3);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsG],flipud(GROIaveMaxFilt));
                        %plot heatmap for times of intensities for each
                        %region
                        title('green maximum intensity');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                        
                        if numROIsG == numROIsR
                            subplot(4,1,4);
                            overlayIm = zeros([size(RROIaveMaxFilt) 3]);
                            overlayIm(:,:,1) = flipud((RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./...
                                (max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt))));
                            overlayIm(:,:,2) = flipud((GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./...
                                (max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt))));
                            image(OffsetRotMatch(:,1),[1:numROIsR],overlayIm);
                            xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                            xlabel('time (sec)');
                        end
                    end
                    colormap(brewermap(64, 'Blues'));

                    set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                    print(act,strcat(ROIfName,'_BumpAct'),'-dpdf'); 

                    delete(act);

                    %clear fullpath stackMaxInt stackReg stackXYfiltMax ROIaveMax positionDat;
                    %Kris no longer clear this, want to have a look at it
                end
            end
        end
    end
    display('finished')
end
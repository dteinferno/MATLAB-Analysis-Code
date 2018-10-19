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

%% Step through each imaging stack in the directory, performing image registration and saving the shift coordinates

% Specify the number of pixels in from the edge to use for image
% registration
Px = 10;

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
    
    % Specify the path where the files will be saved
    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
    if exist(newPath,'dir') == 0
        mkdir(newPath);
    end
    
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif') & ~contains(tifName,'Anatomy')
            
            % Set the filename
            fName = strcat(newPath,tifName(1:end-4));

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
            if allPathname{dirNow}.numColors == 1
                numVolumes = floor(num_images/num_planes);
                
                % Find the max intensity projection
                stackMaxInt = double(zeros(height,width,numVolumes));
                miniStack = double(zeros(height,width,num_planes-num_discards));
                h = waitbar(0.0,'Loading TIFF stack (max calc)...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Loading TIFF stack...');
                for incIm = 1:num_images
                    if mod(incIm,100)==0
                        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
                    end
                    if mod(incIm,num_planes) ~= 0 & mod(incIm,num_planes) < num_planes-num_discards+1
                        miniStack(:,:,mod(incIm,num_planes)) = vol(:,:,incIm)';
                        if mod(incIm,num_planes) == num_planes-num_discards
                            stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
                        end
                    end
                end
                delete(h);
                
                % Register the stacks
                [stackReg, regCoords] = imRegSimple(stackMaxInt, Px); % Image correct the stacks
            else
                numVolumes = floor(num_images/num_planes/2);
                
                % Find the max intensity projection
                RstackMaxInt = double(zeros(height,width,numVolumes));
                stackMaxInt = double(zeros(height,width,numVolumes));
                RminiStack = double(zeros(height,width,(num_planes-num_discards)));
                GminiStack = double(zeros(height,width,(num_planes-num_discards)));
                h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Loading TIFF stack...');
                for incIm = 1:num_images
                    if mod(incIm,100)==0
                        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
                    end
                    if mod(ceil(incIm/2),num_planes) ~= 0
                        if mod(incIm,2)
                            RminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                        else
                            GminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                        end
                        if (mod(incIm/2,num_planes) == num_planes-num_discards)
                            RstackMaxInt(:,:,ceil(incIm/(2*num_planes))) = max(RminiStack,[],3);
                            stackMaxInt(:,:,ceil(incIm/(2*num_planes))) = max(GminiStack,[],3);
                        end
                    end
                end
                delete(h);
                
                % Register the stacks
                [stackReg RstackReg, regCoords] = imRegSimple2Color(stackMaxInt, RstackMaxInt, Px); % Image correct the stacks
            end
            
            % Find the mean registered images (for use specifyng ROIs
            % later)
            meanOrigStack = mean(stackMaxInt,3);
            meanRegStack = mean(stackReg,3);
            if allPathname{dirNow}.numColors == 2
                meanOrigStack_R = mean(RstackMaxInt,3);
                meanRegStack_R = mean(RstackReg,3);
            end
            
            % Save the appropriate data
            if allPathname{dirNow}.numColors == 1
                save(strcat(fName,'_RegDat.mat'), 'Px','regCoords','meanRegStack');
            else
                save(strcat(fName,'_RegDat.mat'), 'Px','regCoords','meanRegStack','meanRegStack_R');
            end
            
            reg = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,2,1);
            imagesc(meanOrigStack)
            axis equal;
            axis off;
            title('stack as taken');

            subplot(2,2,2);
            imagesc(meanRegStack)
            axis equal;
            axis off;
            title('registered stack');
            
            if allPathname{dirNow}.numColors == 2
                subplot(2,2,3);
                imagesc(meanOrigStack_R)
                axis equal;
                axis off;
                title('stack as taken');

                subplot(2,2,4);
                imagesc(meanRegStack_R)
                axis equal;
                axis off;
                title('registered stack');
            end

            colormap(brewermap(64, 'Blues'));
            
            set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
            print(reg,strcat(fName,'_Registered'),'-dpdf'); 
            
            delete(reg);
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
clear stackMaxInt;

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
    
    % Specify the path where the files are saved
    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
    
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif') & ~contains(tifName,'Anatomy')
            tifName
            
            % Get the filename and load the data
            fName = strcat(newPath,tifName(1:end-4));
            load(strcat(fName,'_RegDat.mat'));
            
            stackMaxInt{1} = meanRegStack;
            if exist('meanRegStack_R')
                stackMaxInt{2} = meanRegStack_R;
            end
               
            ROIs = SetROIsClean(allPathnameNow,tifName,allPathname{dirNow}.imRegion,allPathname{dirNow}.numColors, stackMaxInt);
            
            newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
            ROIfName = strcat(newPath,tifName(1:end-4));
    
            save(strcat(ROIfName,'_ROIs','.mat'), 'ROIs');
            clear ROIs;
        end
    end 
end

%% For each stack calculate ROIs

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
                    
                    clear stackMaxInt stackReg;
                    
                    % Get the maximum intensity stack(s)
                    if allPathname{dirNow}.numColors == 1
                        [stackMaxIntNow, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,0);
                        stackMaxInt{1} = stackMaxIntNow;
                    elseif allPathname{dirNow}.numColors == 2
                        [RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] =...
                            ImDatLoadBigtiff2Color(tifName,allPathnameNow,0);
                        stackMaxInt{1} = GstackMaxInt;
                        stackMaxInt{2} = RstackMaxInt;
                    end
                    
                    % Image correct the stacks:
                    % Load the registration coordinates
                    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
                    ROIfName = strcat(newPath,tifName(1:end-4));
                    load(strcat(ROIfName,'_RegDat','.mat'));
                    
                    % Create blank stacks for the registered images
                    imWidth = size(stackMaxInt{1},1);
                    imHeight = size(stackMaxInt{1},2);
                    stackReg{1} = zeros(size(stackMaxInt{1}));
                    if allPathname{dirNow}.numColors == 2
                        stackReg{2} = zeros(size(stackMaxInt{2}));
                    end

                    % Specify the registration image and the area to be registered
                    regArea = [Px Px size(stackMaxInt{1},1)-Px*2+1 size(stackMaxInt{1},2)-Px*2+1];

                    % relative offset of position of subimages
                    rect_offset = [Px
                                   Px];

                    h = waitbar(0.0,'Registering stack...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Registering images...');

                    % Register each image
                    for imNow = 1:size(stackMaxInt{1},3)

                        waitbar(imNow/size(stackMaxInt{1},3),h,['Loading frame# ' num2str(imNow) ' out of ' num2str(size(stackMaxInt{1},3))]);

                        for colNow = 1:allPathname{dirNow}.numColors
                        % Choose subregion of image
                            current_im = stackMaxInt{colNow}(:,:,imNow);
                            sub_EB = imcrop(current_im,regArea);

                            % coordinates for placing it in the registered stack
                            xbegin = regCoords(1,imNow);
                            xend   = regCoords(2,imNow);
                            ybegin = regCoords(3,imNow);
                            yend   = regCoords(4,imNow);

                            % put it in the stack
                            stackMaxIntRegTMP = zeros(imWidth+2*Px,imHeight+2*Px);
                            stackMaxIntRegTMP(ybegin:2*Px+yend,xbegin:2*Px+xend) = stackMaxInt{colNow}(:,:,imNow);
                            stackReg{colNow}(:,:,imNow) = imcrop(stackMaxIntRegTMP,[Px+1 Px+1 imWidth-1 imHeight-1]);
                        end

                    end

                    delete(h);
                                         
                    % Rotate the stacks
                    if strcmp(allPathname{dirNow}.imRegion{1},'EB') || strcmp(allPathname{dirNow}.imRegion{1},'FB')
                        rotAng = -70;
                    else
                        rotAng = 90;
                    end
                    
                    if allPathname{dirNow}.numColors == 1
                        stkMean = mean(stackReg{1},3);
                        stkMean = imrotate(stkMean,rotAng);
                        stackMaxIntRot = zeros([size(stkMean) size(stackReg{1},3)]);
                        for frame = 1:size(stackReg{1},3)
                            stackMaxIntRot(:,:,frame) = imrotate(stackReg{1}(:,:,frame),rotAng);
                        end
                    else
                        for colNow = 1:2
                            stkMean = mean(stackReg{colNow},3);
                            stkMean = imrotate(stkMean,rotAng);
                            stackMaxIntRot = zeros([size(stkMean) size(stackReg{colNow},3)]);
                            for frame = 1:size(stackReg{colNow},3)
                                stackMaxIntRot(:,:,frame) = imrotate(stackReg{colNow}(:,:,frame),rotAng);
                            end
                            if colNow == 1
                                GstackMaxIntRot = stackMaxIntRot;
                            else
                                RstackMaxIntRot = stackMaxIntRot;
                            end
                        end
                    end

                    % Load the movement data
                    positionDat = VRDatLoad(moveName,allPathnameNow,0);
                    
                    % Load the ROIs
                    load(strcat(ROIfName,'_ROIs','.mat'));
                    numROIsG=size(ROIs{1},1);
                    if allPathname{dirNow}.numColors == 2
                        numROIsR=size(ROIs{2},1);
                    end
                    if allPathname{dirNow}.numColors == 1
                        numPlanes = length(positionDat.tFrameGrab)/length(stackMaxInt{1});
                    elseif allPathname{dirNow}.numColors == 2
                        numPlanes = length(positionDat.tFrameGrab)/length(RstackMaxInt);
                    end

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
                        RROIs = squeeze(ROIs{2});
                        GROIs = squeeze(ROIs{1});
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
                    
                    
                    if allPathname{dirNow}.numColors == 1
                        subplot(4,1,2);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsG],flipud(ROIaveMaxFilt));
                        title('maximum intensity');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                        
                        vRot = diff(pi/180*PosRotMatch)./mean(diff(OffsetRotMatch(:,1)));
                        subplot(4,1,3);
                        plot(OffsetRotMatch(2:end,1),vRot,'b');
                        title('vRot');
                        ylim([-pi pi]);
                        ylabel('vRot (rad/s)');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                        
                        vF = sqrt(diff(PosForMatch).^2+diff(PosLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));
                        subplot(4,1,4);
                        plot(OffsetRotMatch(2:end,1),vF,'Color',[1 0.5 0]);
                        title('vF');
                        xlabel('time (sec)');
                        ylim([0 2]);
                        ylabel('vF (cm/s)');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                        
                    elseif allPathname{dirNow}.numColors == 2
                        subplot(4,1,2);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsR],flipud(RROIaveMaxFilt));
                        title('red maximum intensity');
                        xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
               
                        subplot(4,1,3);
                        imagesc(OffsetRotMatch(:,1),[1:numROIsG],flipud(GROIaveMaxFilt));
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

                    clear fullpath stackMaxInt stackReg stackXYfiltMax ROIaveMax positionDat;
                end
            end
        end
    end    
end
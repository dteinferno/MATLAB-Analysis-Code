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
            
            % Find the mean registered images (for use specifyng ROIs
            % later)
            meanOrigStack = mean(stackMaxInt,3);
            meanRegStack = mean(stackReg,3);
            
            % Save the appropriate data
            save(strcat(fName,'_RegDat.mat'), 'Px','regCoords','meanRegStack');
            
            reg = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(1,2,1);
            imagesc(meanOrigStack)
            axis equal;
            axis off;
            title('stack as taken');

            subplot(1,2,2);
            imagesc(meanRegStack)
            axis equal;
            axis off;
            title('registered stack');

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

% Specify the rotation angle
rotAng = -70;

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
            
            % Get the filename and load the data
            fName = strcat(newPath,tifName(1:end-4));
            load(strcat(fName,'_RegDat.mat'));
            
            tifName

            % Load the reference ROIs
            newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
            pathParts = strsplit(tifName,'_');
            IDLoc = strfind(tifName,pathParts{end-1});
            ROIfName = strcat(newPath,tifName(1:IDLoc-2));
            load(strcat(ROIfName,'_ReferenceROIs'));

            numROIsG = 16;
            meanRegRot = imrotate(meanRegStack,rotAng);

            stackSize = size(meanRegRot);

            ROIs{1} = zeros([numROIsG, stackSize]);

            % Show the stack
            A = meanRegRot;
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);

            % Find the change in angles
            deltaAng = 2*pi/numROIsG;

            % Specify an outer ellipse
            hOut = impoly(gca,position{1});
            positionOut = wait(hOut);
            ellipseCenterOut = [mean(positionOut(:,1)) mean(positionOut(:,2))];
            ellipseHeightOut = max(positionOut(:,2)) - min(positionOut(:,2));
            ellipseWidthOut = max(positionOut(:,1)) - min(positionOut(:,1));

            % Specify an inner ellipse
            hIn = impoly(gca,position{2});
            positionIn = wait(hIn);
            ellipseCenterIn = [mean(positionIn(:,1)) mean(positionIn(:,2))];
            ellipseHeightIn = max(positionIn(:,2)) - min(positionIn(:,2));
            ellipseWidthIn = max(positionIn(:,1)) - min(positionIn(:,1));

            for incROI = 1:numROIsG
                for angs=1:5
                    xROI(angs) = ellipseCenterOut(1)+ellipseWidthOut/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                    yROI(angs) = ellipseCenterOut(2)+ellipseHeightOut/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                end
                for angs=6:10
                    xROI(16-angs) = ellipseCenterIn(1)+ellipseWidthIn/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
                    yROI(16-angs) = ellipseCenterIn(2)+ellipseHeightIn/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
                end
                ROIs{1}(incROI,:,:) = roipoly(A, xROI,yROI);
            end
            delete(hf);
            save(strcat(fName,'_ROIs','.mat'), 'ROIs');
            clear ROIs;
        end
    end 
end

%% For each stack perform the following operations: register, calculate ROIs

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
                    
                    % Get the maximum intensity stack
                    [stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,0);
                    
                    % Image correct the stacks:
                    % Load the registration coordinates
                    newPath = strrep(allPathnameNow,'Downloads\Data','Documents\RawAnalysis');
                    ROIfName = strcat(newPath,tifName(1:end-4));
                    load(strcat(ROIfName,'_RegDat','.mat'));
                    
                    % Create blank stacks for the registered images
                    imWidth = size(stackMaxInt,1);
                    imHeight = size(stackMaxInt,2);
                    stackReg = zeros(size(stackMaxInt));

                    % Specify the registration image and the area to be registered
                    ref_im = mean(stackMaxInt,3);
                    orig_res = [0 0 size(stackMaxInt,1) size(stackMaxInt,2)];
                    regArea = [Px Px size(stackMaxInt,1)-Px*2+1 size(stackMaxInt,2)-Px*2+1];

                    % relative offset of position of subimages
                    rect_offset = [Px
                                   Px];

                    h = waitbar(0.0,'Registering stack...');
                    set(h,'Position',[50 50 360 72]);
                    set(h,'Name','Registering images...');

                    % Register each image
                    for imNow = 1:size(stackMaxInt,3)

                        waitbar(imNow/size(stackMaxInt,3),h,['Loading frame# ' num2str(imNow) ' out of ' num2str(size(stackMaxInt,3))]);

                        % Choose subregion of image
                        current_im = stackMaxInt(:,:,imNow);
                        sub_EB = imcrop(current_im,regArea);

                        % coordinates for placing it in the registered stack
                        xbegin = regCoords(1,imNow);
                        xend   = regCoords(2,imNow);
                        ybegin = regCoords(3,imNow);
                        yend   = regCoords(4,imNow);

                        % put it in the stack
                        stackMaxIntRegTMP = zeros(imWidth+2*Px,imHeight+2*Px);
                        stackMaxIntRegTMP(ybegin:2*Px+yend,xbegin:2*Px+xend) = stackMaxInt(:,:,imNow);
                        stackReg(:,:,imNow) = imcrop(stackMaxIntRegTMP,[Px+1 Px+1 imWidth-1 imHeight-1]);

                    end

                    delete(h);
                                         
                    % Rotate the stacks
                    if strcmp(allPathname{dirNow}.imRegion{1},'EB') || strcmp(allPathname{dirNow}.imRegion{1},'FB')
                        rotAng = -70;
                    else
                        rotAng = 90;
                    end
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
%% Get an ROI from just the first few volumes of a tif
function [ROIs] = SetROIs(allPathname,tifName,imRegion,numColors)

     % Get the tiff header info
    fullpath = strcat(allPathname,tifName);
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
    if numColors == 2
        numVolumes = numVolumes/2;
    end
    
    % Find the max intensity projection
    if numColors == 1
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
                miniStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                if (mod(incIm,num_planes) == num_planes-num_discards)
                    stackMaxInt(:,:,ceil(incIm/(num_planes))) = max(miniStack,[],3);
                end
            end
        end
        delete(h);
    elseif numColors == 2
        RstackMaxInt = double(zeros(height,width,numVolumes));
        GstackMaxInt = double(zeros(height,width,numVolumes));
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
                    GstackMaxInt(:,:,ceil(incIm/(2*num_planes))) = max(GminiStack,[],3);
                end
            end
        end
        delete(h);
    end
    
     % Rotate the stack if the EB was imaged
    if strcmp(imRegion{1},'EB') || strcmp(imRegion{1},'FB')
        rotAng = -70;
    else
        rotAng = 90;
    end
    if numColors == 1
        % Rotate the stacks so that the EB is aligned
        stackMean = mean(stackMaxInt,3);
        stackMean = imrotate(stackMean,rotAng);
        stackMaxIntRot = zeros([size(stackMean) size(stackMaxInt,3)]);
        for frame = 1:size(stackMaxIntRot,3)
            stackMaxIntRot(:,:,frame) = imrotate(stackMaxInt(:,:,frame),rotAng);
        end
    elseif numColors == 2
        RMean = mean(RstackMaxInt,3);
        RMean = imrotate(RMean,rotAng);
        GMean = mean(GstackMaxInt,3);
        GMean = imrotate(GMean,rotAng);
        RstackMaxIntRot = zeros([size(RMean) size(RstackMaxInt,3)]);
        GstackMaxIntRot = zeros([size(GMean) size(GstackMaxInt,3)]);
        for frame = 1:size(RstackMaxIntRot,3)
            RstackMaxIntRot(:,:,frame) = imrotate(RstackMaxInt(:,:,frame),rotAng);
            GstackMaxIntRot(:,:,frame) = imrotate(GstackMaxInt(:,:,frame),rotAng);
        end
    end

    % Load the reference ROIs
    newPath = strrep(allPathname,'Downloads\Data','Documents\RawAnalysis');
    pathParts = strsplit(tifName,'_');
    IDLoc = strfind(tifName,pathParts{end-1});
    ROIfName = strcat(newPath,tifName(1:IDLoc-2));
    load(strcat(ROIfName,'_ReferenceROIs'));
    
    % Get the ROIs
    if strcmp(imRegion{1},'EB')
        numROIsG = 16;
    elseif strcmp(imRegion{1},'PB')
        numROIsG = 18;
    elseif strcmp(imRegion{1},'FB')
        numROIsG = 8;
    elseif strcmp(imRegion{1},'NO')
        numROIsG = 2;
    elseif strcmp(imRegion{1},'other')
        if numColors == 1
            numROIsG = length(position);
        else
            numROIsG = length(positionG);
        end
    elseif 1
        disp('Functionality not yet enabled');
    end
    if numColors == 2
        if strcmp(imRegion{2},'EB')
            numROIsR = 16;
        elseif strcmp(imRegion{2},'PB')
            numROIsR = 18;
        elseif strcmp(imRegion{2},'FB')
            numROIsR = 8;
        elseif strcmp(imRegion{2},'NO')
            numROIsR = 2;
        elseif strcmp(imRegion{2},'other')
            numROIsR = length(positionR);
        end
    end
    if numColors == 1
        stackSize = [size(stackMaxIntRot,1) size(stackMaxIntRot,2)];
    else
        stackSize = [size(RstackMaxIntRot,1) size(RstackMaxIntRot,2)];
    end
    ROIs{1} = zeros([numROIsG, stackSize]);
    if numColors == 2
        ROIs{2} = zeros([numROIsR, stackSize]);
    end
    
    if numColors == 1
        % Show the stack
        A = mean(stackMaxIntRot,3);
        hf = figure;
        hold;
        imshow(A,[0 max(max(A))]);
    
        % Get the ROI
        if strcmp(imRegion{1},'EB')
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
        else
            for pgs = 1:numROIsG
                pgons{pgs} = impoly(gca,position{pgs});
            end
            wait(pgons{1});
            for pgs = 1:numROIsG
                posNow = getPosition(pgons{pgs});
                ROIs{1}(pgs,:,:) = roipoly(A,posNow(:,1),posNow(:,2));
            end
        end
        delete(hf);
    elseif numColors == 2
        if strcmp(imRegion{1},'EB')
            % Find the change in angles
            deltaAng = 2*pi/numROIsG;
            
             % Get the ROIs for the green channel
            A = mean(GstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            % Specify an outer ellipse
            hOut = impoly(gca,positionG{1});
            positionOut = wait(hOut);
            ellipseCenterOut = [mean(positionOut(:,1)) mean(positionOut(:,2))];
            ellipseHeightOut = max(positionOut(:,2)) - min(positionOut(:,2));
            ellipseWidthOut = max(positionOut(:,1)) - min(positionOut(:,1));

            % Specify an inner ellipse
            hIn = impoly(gca,positionG{2});
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
            
        elseif strcmp(imRegion{1},'PB')
            A = mean(RstackMaxIntRot,3);
            Rim = zeros(size(A,1),size(A,2),3);
            Rim(:,:,1) = (A-min(min(A)))./(max(max(A))-min(min(A)));
            B = mean(GstackMaxIntRot,3);
            Gim = zeros(size(B,1),size(B,2),3);
            Gim(:,:,2) = (B-min(min(B)))./(max(max(B))-min(min(B)));
            overlayIm = Rim+Gim;
            
            hf = figure;
            subplot(2,3,1);
            imshow(Rim);
            subplot(2,3,4);
            imshow(Gim);
            subplot(2,3,[2:3 5:6]);
            hold;
            imshow(overlayIm);
            for pgs = 1:numROIsR
                pgons{pgs} = impoly(gca,positionR{pgs});
            end
            wait(pgons{1});
            for pgs = 1:numROIsG
                position = getPosition(pgons{pgs});
                ROIs{1}(pgs,:,:) = roipoly(2*overlayIm,position(:,1),position(:,2));
                ROIs{2}(pgs,:,:) = roipoly(2*overlayIm,position(:,1),position(:,2));
            end
            delete(hf);
            
        elseif strcmp(imRegion{1},'FB') || strcmp(imRegion{1},'other')
            A = mean(GstackMaxIntRot,3);
            
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            for pgs = 1:numROIsG
                pgons{pgs} = impoly(gca,positionG{pgs});
            end
            wait(pgons{1});
            for pgs = 1:numROIsG
                position = getPosition(pgons{pgs});
                ROIs{1}(pgs,:,:) = roipoly(A,position(:,1),position(:,2));
            end
            delete(hf);
            
        elseif 1
            disp('Functionality not yet enabled');
        end
        
        if strcmp(imRegion{2},'EB')
            % Find the change in angles
            deltaAng = 2*pi/numROIsR;
            
            % Show the red channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
        
            % Specify an outer ellipse
            hOut = impoly(gca,positionR{1});
            positionOut = wait(hOut);
            ellipseCenterOut = [mean(positionOut(:,1)) mean(positionOut(:,2))];
            ellipseHeightOut = max(positionOut(:,2)) - min(positionOut(:,2));
            ellipseWidthOut = max(positionOut(:,1)) - min(positionOut(:,1));

            % Specify an inner ellipse
            hIn = impoly(gca,positionR{2});
            positionIn = wait(hIn);
            ellipseCenterIn = [mean(positionIn(:,1)) mean(positionIn(:,2))];
            ellipseHeightIn = max(positionIn(:,2)) - min(positionIn(:,2));
            ellipseWidthIn = max(positionIn(:,1)) - min(positionIn(:,1));
            
            for incROI = 1:numROIsR
                for angs=1:5
                    xROI(angs) = ellipseCenterOut(1)+ellipseWidthOut/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                    yROI(angs) = ellipseCenterOut(2)+ellipseHeightOut/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                end
                for angs=6:10
                    xROI(16-angs) = ellipseCenterIn(1)+ellipseWidthIn/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
                    yROI(16-angs) = ellipseCenterIn(2)+ellipseHeightIn/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
                end
                ROIs{2}(incROI,:,:) = roipoly(A, xROI,yROI);
            end
            
            delete(hf);
        elseif strcmp(imRegion{2},'FB') || strcmp(imRegion{2},'other')
            % Show the red channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            for pgs = 1:numROIsR
                pgons{pgs} = impoly(gca,positionR{pgs});
            end
            wait(pgons{1});
            for pgs = 1:numROIsR
                position = getPosition(pgons{pgs});
                ROIs{2}(pgs,:,:) = roipoly(A,position(:,1),position(:,2));
            end
            delete(hf);
        elseif strcmp(imRegion{2},'NO')
            A = mean(RstackMaxIntRot,3);
            
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            for pgs = 1:numROIsR
                pgons{pgs} = impoly(gca,positionR{pgs});
            end
            wait(pgons{1});
            for pgs = 1:numROIsR
                position = getPosition(pgons{pgs});
                ROIs{2}(pgs,:,:) = roipoly(A,position(:,1),position(:,2));
            end
            delete(hf);
        end
    end
end
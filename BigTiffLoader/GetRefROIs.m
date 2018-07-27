%% Get ROIs for the E-PGs and ring neurons
function GetRefROIs(imRegion,numColors)

    % Get the 2P info
    [tifName,allPathname] = uigetfile('*.tif','Select the Two Photon Data');
    fullpath = strcat(allPathname,tifName);
    cd(allPathname);

    % Get the tiff header info

    fullpath = strcat(allPathname,tifName);
    reader = ScanImageTiffReader(fullpath);
    desc = reader.metadata;
    planeLoc = strfind(desc, 'numFramesPerVolume');
    discardLoc = strfind(desc, 'numDiscardFlybackFrames');
    num_planes = str2num(desc(planeLoc+21:planeLoc+22));
    num_discards = str2num(desc(discardLoc+26));

    % Load the tifStack
    fullpath;
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

    % Rotate the stack if the EB or FB was imaged
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

    % Show the average MIP and get ROIs
    if numColors == 1
        position = {};
        
        % Show the stack
        A = mean(stackMaxIntRot,3);
        hf = figure;
        hold;
        imshow(A,[0 max(max(A))]);
    
        % Get the ROI
        if strcmp(imRegion{1},'EB')
            % Specify an outer ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            position{1} = posNow;

            % Specify an inner ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            position{2} = posNow;
        elseif strcmp(imRegion{1},'PB')
            for pgs = 1:18
                h = impoly;
                posNow = wait(h);
                position{pgs} = posNow;
            end
        elseif strcmp(imRegion{1},'FB')
            for pgs = 1:8
                h = impoly;
                posNow = wait(h);
                position{pgs} = posNow;
            end
        elseif strcmp(imRegion{1},'NO')
            % Specify a left ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            position{1} = posNow;

            % Specify a right ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            position{2} = posNow;
        elseif strcmp(imRegion{1},'other')
            num_ellipses = input('Number of ellipses?');
            num_pgons = input('Number of polygons?');
            num_ROIs = num_ellipses+num_pgons;

            ROIs = zeros(num_ROIs,size(stackMaxInt,1),size(stackMaxInt,2));
            ROIverts = cell(num_ROIs,1);
            for els = 1:num_ellipses
                h = imellipse;
                posNow = wait(h);
                position{els} = posNow;
            end
            for pgs = 1:num_pgons
                h = impoly;
                posNow = wait(h);
                position{num_ellipses+pgs} = posNow;
            end
        elseif 1
            disp('Functionality not yet enabled');
        end
        delete(hf);
    elseif numColors == 2
        
        positionG = {};
        positionR = {};
        if strcmp(imRegion{1},'EB')
            % Get the ROIs for the green channel
            A = mean(GstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            % Specify an outer ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            positionG{1} = posNow;

            % Specify an inner ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            positionG{2} = posNow;
            
            delete(hf);
        elseif strcmp(imRegion{1},'PB')
            num_glomeruli = 18;
            pgons = {};
            scaleFact = round(size(GstackMaxIntRot,1)/70);
            defaultGlomShape = [1,0;3,0;4,1;4,4;3,5;1,5;0,4;0,1];
            defaultGlomShape = scaleFact*defaultGlomShape;
            glomOffsets = [-7.5,2;-7.5,0.5;-6.5,-0.5;-5.5,-1.5;-4.5,-1.5;-3.5,-1.25;-2.5,-1;-1.5,-0.75;-0.5,-0.5];
            glomOffsets = vertcat(glomOffsets,flipud(glomOffsets.*[-1,1]));
            glomOffsets = 0.5*[size(GstackMaxIntRot,1),size(GstackMaxIntRot,2)]+...
                4*scaleFact*glomOffsets;
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
            imshow(overlayIm);
            
            RMult = input('Red channel multiplier? ');
            GMult = input('Green channel multiplier? ');
            
            subplot(2,3,1);
            cla;
            imshow(RMult*Rim);
            subplot(2,3,4);
            cla;
            imshow(GMult*Gim);
            subplot(2,3,[2:3 5:6]);
            cla;
            hold;
            imshow(RMult*Rim+GMult*Gim);
            
            for pgs = 1:num_glomeruli
                pgons{pgs} = impoly(gca,defaultGlomShape+glomOffsets(pgs,:));
                api = iptgetapi(pgons{pgs});
                api.setColor([(pgs-1)/(num_glomeruli-1) 1 1-(pgs-1)/(num_glomeruli-1)]);
            end
            wait(pgons{1});
            
            for pgs = 1:num_glomeruli
                position = getPosition(pgons{pgs});
                if contains(tifName,'Anatomy')
                    position = 256./size(RstackMaxIntRot,1).*position;
                end
                
                positionR{pgs} = position;
                positionG{pgs} = position;
            end      
            delete(hf);
        elseif strcmp(imRegion{1},'FB')
            % Get the ROIs for the green channel
            A = mean(GstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            for pgs = 1:8
                h = impoly;
                posNow = wait(h);
                positionG{pgs} = posNow;
            end
            delete(hf);
        elseif strcmp(imRegion{1},'other')
            % Get the ROIs for the green channel
            A = mean(GstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            num_ellipses = input('Number of ellipses?');
            num_pgons = input('Number of polygons?');

            for els = 1:num_ellipses
                h = imellipse;
                posNow = wait(h);
                positionG{els} = posNow;
            end
            for pgs = 1:num_pgons
                h = impoly;
                posNow = wait(h);
                positionG{num_ellipses+pgs} = posNow;
            end
            delete(hf);
        end
            
        if strcmp(imRegion{2},'EB')
            % Show the red channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
        
            % Specify an outer ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            positionR{1} = posNow;

            % Specify an inner ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            positionR{2} = posNow;
            
            delete(hf);
        elseif strcmp(imRegion{2},'FB')
            % Get the ROIs for the green channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            for pgs = 1:8
                h = impoly;
                posNow = wait(h);
                positionR{pgs} = posNow;
            end
            delete(hf);
        elseif strcmp(imRegion{2},'NO')
            % Get the ROIs for the green channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            % Specify a left ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            positionR{1} = posNow;

            % Specify a right ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            positionR{2} = posNow;
            
            delete(hf);
        elseif strcmp(imRegion{2},'other')
            % Get the ROIs for the red channel
            A = mean(RstackMaxIntRot,3);
            hf = figure;
            hold;
            imshow(A,[0 max(max(A))]);
            
            num_ellipses = input('Number of ellipses?');
            num_pgons = input('Number of polygons?');

            for els = 1:num_ellipses
                h = imellipse;
                posNow = wait(h);
                positionR{els} = posNow;
            end
            for pgs = 1:num_pgons
                h = impoly;
                posNow = wait(h);
                positionR{num_ellipses+pgs} = posNow;
            end
        end
    end
    
    newPath = strrep(allPathname,'Downloads\Data','Documents\RawAnalysis');
    if exist(newPath,'dir') == 0
        mkdir(newPath);
    end
    pathParts = strsplit(tifName,'_');
    IDLoc = strfind(tifName,pathParts{end-1});
    ROIfName = strcat(newPath,tifName(1:IDLoc-2));
    if numColors == 1
        save(strcat(ROIfName,'_ReferenceROIs'),'position');
    else
        save(strcat(ROIfName,'_ReferenceROIs'),'positionR','positionG');
    end
end
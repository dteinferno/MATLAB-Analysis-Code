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
        
        % Show the stack
        A = mean(stackMaxIntRot,3);
        hf = figure('units','normalized','outerposition',[0 0 1 1]);
        hold;
        imshow(A,[0 max(max(A))]);
        axis equal;
        axis off;
    
        position = OneColorROIs(stackMaxIntRot, imRegion, {});
        
        done = input('Fix the ROIs (y or n)');
        
        while strcmp(done,'y')
            cla;
            imshow(A,[0 max(max(A))]);
            axis equal;
            axis off;
            
            position = OneColorROIs(stackMaxIntRot, imRegion, position);
            
            done = input('Fix the ROIs (y or n)');
        end
        
        delete(hf);
    elseif numColors == 2
        [positionG, positionR] = TwoColorROIs(GstackMaxIntRot, RstackMaxIntRot,imRegion);
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
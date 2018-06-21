%% Get an ROI from just the first few volumes of a tif
function [RROIs,GROIs] = Get2ColorROIs(allPathname,tifName,rotAng)

% Get the tiff header info
    fullpath = strcat(allPathname,tifName);
    reader=ScanImageTiffReader(fullpath);
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
    numVolumes = floor(num_images/num_planes/2);

    % Find the max intensity projection
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
    
    % Rotate the stacks so that the EB is aligned
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

    % Get the ROIs for the red channel
    A = mean(RstackMaxIntRot,3);
    hf = figure;
    hold;
    imshow(A,[0 max(max(A))]);
    h = imellipse;
    position = wait(h);
%     position(length(position)+1,:) = position(1,:);
%     ellipseCenter = [position(7,1) position(1,2)];
%     ellipseHeight = position(7,2) - position(19,2);
%     ellipseWidth = position(1,1) - position(13,1);
    ellipseCenter = [mean(position(:,1)) mean(position(:,2))];
    ellipseHeight = max(position(:,2)) - min(position(:,2));
    ellipseWidth = max(position(:,1)) - min(position(:,1));
    num_ROIs = 16;
    deltaAng = 2*pi/num_ROIs;
    delete(hf);

    RROIs = zeros(num_ROIs,size(RstackMaxIntRot,1),size(RstackMaxIntRot,2));
    for incROI = 1:num_ROIs
        xROI(1) = ellipseCenter(1);
        yROI(1) = ellipseCenter(2);
        for angs=1:5
            xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        end
        RROIs(incROI,:,:) = roipoly(A, xROI,yROI);
    end

    % Get the ROIs for the green channel
    A = mean(GstackMaxIntRot,3);
    hf = figure;
    hold;
    imshow(A,[0 max(max(A))]);
    
    num_ellipses = input('Number of ellipses?');
    num_pgons = input('Number of polygons?');
    num_ROIs = num_ellipses+num_pgons;

    GROIs = zeros(num_ROIs,size(GstackMaxIntRot,1),size(GstackMaxIntRot,2));
    for els = 1:num_ellipses
        h = imellipse;
        position = wait(h);
        GROIs(els,:,:) = roipoly(A,position(:,1),position(:,2));
    end
    for pgs = 1:num_pgons
        h = impoly;
        position = wait(h);
        GROIs(num_ellipses+pgs,:,:) = roipoly(A,position(:,1),position(:,2));
    end
    
    delete(hf);
    
%     h = imellipse;
%     position = wait(h);
% %     position(length(position)+1,:) = position(1,:);
% %     ellipseCenter = [position(7,1) position(1,2)];
% %     ellipseHeight = position(7,2) - position(19,2);
% %     ellipseWidth = position(1,1) - position(13,1);
%     ellipseCenter = [mean(position(:,1)) mean(position(:,2))];
%     ellipseHeight = max(position(:,2)) - min(position(:,2));
%     ellipseWidth = max(position(:,1)) - min(position(:,1));
%     num_ROIs = 16;
%     deltaAng = 2*pi/num_ROIs;
%     delete(hf);
% 
%     GROIs = zeros(num_ROIs,256,256);
%     for incROI = 1:num_ROIs
%         xROI(1) = ellipseCenter(1);
%         yROI(1) = ellipseCenter(2);
%         for angs=1:5
%             xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
%             yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
%         end
%         GROIs(incROI,:,:) = roipoly(A, xROI,yROI);
%     end
end
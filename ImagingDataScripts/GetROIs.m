%% Get an ROI from just the first few volumes of a tif
function [ROIs] = GetROIs(allPathname,tifName)

    % Specify the number of volumes to look at
%     numVolumes = 100;

    % Get the 2P info;
    fullpath = strcat(allPathname,tifName);
    info = imfinfo(fullpath);
    recSpecs = info(1).ImageDescription;
    planeLoc = strfind(recSpecs, 'numFramesPerVolume');
    discardLoc = strfind(recSpecs, 'numDiscardFlybackFrames');
    num_planes = str2num(recSpecs(planeLoc+21));
    num_discards = str2num(recSpecs(discardLoc+26));
    if isempty(num_planes)
        hTiff = Tiff(fullpath);
        configStr = hTiff.getTag('Software');
        planeLoc = strfind(configStr, 'numFramesPerVolume');
        discardLoc = strfind(configStr, 'numDiscardFlybackFrames');
        num_planes = str2num(configStr(planeLoc+21:planeLoc+22));
        num_discards = str2num(configStr(discardLoc+26));
    end
%     width = info(1).Width;
%     height = info(1).Height;
    vol=ScanImageTiffReader(fullpath).data();
    width = size(vol,1);
    height = size(vol,2);
    num_images = size(vol,3);

    % Find the max intensity projection
%     stackMaxInt = double(zeros(height,width,numVolumes));
    stackMaxInt = double(zeros(height,width,num_images));
    miniStack = double(zeros(height,width,num_planes-num_discards));
%     for incIm = 1:numVolumes*num_planes
    for incIm = 1:num_images
        if mod(incIm,num_planes) ~= 0
            miniStack(:,:,1+mod(incIm-1,num_planes-num_discards)) = double(imread(fullpath, incIm, 'Info', info));
            if (1+mod(incIm-1,num_planes-num_discards) == num_planes-num_discards)
                stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
            end
        end
    end

    % Get the ROIs
    A = mean(stackMaxInt,3);
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

    ROIs = zeros(num_ROIs,256,256);
    for incROI = 1:num_ROIs
        xROI(1) = ellipseCenter(1);
        yROI(1) = ellipseCenter(2);
        for angs=1:5
            xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        end
        ROIs(incROI,:,:) = roipoly(A, xROI,yROI);
    end
end
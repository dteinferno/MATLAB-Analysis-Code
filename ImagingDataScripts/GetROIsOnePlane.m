%% Get an ROI from just the first few volumes of a tif
function [ROIs] = GetROIsOnePlane(allPathname,tifName, planeNum)

    % Specify the number of volumes to look at
    numVolumes = 100;

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
    width = info(1).Width;
    height = info(1).Height;

    % Get out the designated plane
    Rstack = double(zeros(height,width,numVolumes));
    Gstack = double(zeros(height,width,numVolumes));
    for incIm = 1:numVolumes*num_planes*2
        if mod(2*num_planes+(incIm-2*planeNum),2*num_planes) == 0
            Rstack(:,:,ceil(incIm/(2*num_planes))) = double(imread(fullpath, incIm-1, 'Info', info));
            Gstack(:,:,ceil(incIm/(2*num_planes))) = double(imread(fullpath, incIm, 'Info', info));
        end
    end

    % Get the ROIs
    num_ellipses = 2;%input('Number of ellipses?');
    num_pgons = 0;%input('Number of polygons?');
    num_ROIs = num_ellipses+num_pgons;
    ROIs = zeros(num_ROIs,width,height);
    A = mean(Rstack,3);
    Rim = zeros(size(A,1),size(A,2),3);
    Rim(:,:,1) = (A-min(min(A)))./(max(max(A))-min(min(A)));
    B = mean(Gstack,3);
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
    for els = 1:num_ellipses
        h = imellipse;
        position = wait(h);
        ROIs(els,:,:) = roipoly(A,position(:,1),position(:,2));
    end
    for pgs = 1:num_pgons
        h = impoly;
        position = wait(h);
        ROIs(num_ellipses+pgs,:,:) = roipoly(A,position(:,1),position(:,2));
    end
    delete(hf);
end
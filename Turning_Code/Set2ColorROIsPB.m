%% Get an ROI from just the first few volumes of a tif
function [ROIs] = Set2ColorROIsPB(allPathname,tifName)

    % Specify the number of volumes to look at
    numVolumes = 100;

    % Get the 2P info;
    fullpath = strcat(allPathname,tifName);
    info = imfinfo(fullpath);
    recSpecs = info(1).ImageDescription;
    planeLoc = strfind(recSpecs, 'numFramesPerVolume');
    discardLoc = strfind(recSpecs, 'numDiscardFlybackFrames');
    num_planes = str2num(recSpecs(planeLoc+21:planeLoc+22));
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

    % Find the mean volume projection
    RstackMean = double(zeros(height,width,numVolumes));
    GstackMean = double(zeros(height,width,numVolumes));
    for incIm = 1:numVolumes*num_planes*2
        if mod(ceil(incIm/2),num_planes) ~= 0
            if mod(incIm,2)
                RstackMean(:,:,ceil(incIm/(2*num_planes))) = double(imread(fullpath, incIm, 'Info', info))+ RstackMean(:,:,ceil(incIm/(2*num_planes)));
            else
                GstackMean(:,:,ceil(incIm/(2*num_planes))) = double(imread(fullpath, incIm, 'Info', info))+ GstackMean(:,:,ceil(incIm/(2*num_planes)));
            end
        end
    end

    % Get the ROIs
    num_glomeruli = 18;
	ROIs = zeros(num_glomeruli, width, height);
    A = mean(RstackMean,3);
    Rim = zeros(size(A,1),size(A,2),3);
    Rim(:,:,1) = (A-min(min(A)))./(max(max(A))-min(min(A)));
    B = mean(GstackMean,3);
    Gim = zeros(size(B,1),size(B,2),3);
    Gim(:,:,2) = (B-min(min(B)))./(max(max(B))-min(min(B)));
    overlayIm = Rim+Gim;
    
    load(strcat(tifName(1:5),'ReferenceROIs'));
    
    hf = figure;
    pgons = {};
    subplot(2,3,1);
    imshow(Rim);
    subplot(2,3,4);
    imshow(Gim);
    subplot(2,3,[2:3 5:6]);
    hold;
    imshow(overlayIm);
    for pgs = 1:num_glomeruli
        pgons{pgs} = impoly(gca,pgonPos{pgs});
    end
    wait(pgons{pgs});
    for pgs = 1:num_glomeruli
        position = getPosition(pgons{pgs});
        ROIs(pgs,:,:) = roipoly(2*overlayIm,position(:,1),position(:,2));
    end      
    delete(hf);

end
%% Get an ROI from just the first few volumes of a tif
function [ROIs] = SetROIsPB(allPathname,tifName)

    % Get the 2P info
    fullpath = strcat(allPathname,tifName);
    cd(allPathname);

    % Specify the number of volumes to look at
    numVolumes = 100;

    % Get the 2P info;
    info = imfinfo(fullpath);
    hTiff = Tiff(fullpath);
    configStr = hTiff.getTag('Software');
    planeLoc = strfind(configStr, 'numFramesPerVolume');
    discardLoc = strfind(configStr, 'numDiscardFlybackFrames');
    num_planes = str2num(configStr(planeLoc+21:planeLoc+22));
    num_discards = str2num(configStr(discardLoc+26));

    vol=ScanImageTiffReader(fullpath).data();
    width = size(vol,1);
    height = size(vol,2);
    num_images = size(vol,3);

    % Find the max intensity projection
    stackMaxInt = double(zeros(height,width,num_images));
    miniStack = double(zeros(height,width,num_planes-num_discards));
    for incIm = 1:num_images
        if mod(incIm,num_planes) ~= 0
            miniStack(:,:,1+mod(incIm-1,num_planes-num_discards)) = double(imread(fullpath, incIm, 'Info', info));
            if (1+mod(incIm-1,num_planes-num_discards) == num_planes-num_discards)
                stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
            end
        end
    end

    % Get the ROIs
    num_glomeruli = 18;
	ROIs = zeros(num_glomeruli, width, height);
    A = mean(stackMaxInt,3);
    Rim = zeros(size(A,1),size(A,2),3);
    Rim(:,:,2) = (A-min(min(A)))./(max(max(A))-min(min(A)));
    
    load(strcat(tifName(1:5),'ReferenceROIs'));
    
    hf = figure;
    pgons = {};
    imshow(Rim);
    hold;
    for pgs = 1:num_glomeruli
        pgons{pgs} = impoly(gca,pgonPos{pgs});
    end
    wait(pgons{pgs});
    for pgs = 1:num_glomeruli
        position = getPosition(pgons{pgs});
        ROIs(pgs,:,:) = roipoly(Rim,position(:,1),position(:,2));
    end      
    delete(hf);

end
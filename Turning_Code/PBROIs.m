function PBROIs(numColors)

% Get the 2P info
[imageFilename,imagePathname] = uigetfile('*.tif','Select the Two Photon Data');
fullpath = strcat(imagePathname,imageFilename);
cd(imagePathname);

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

if numColors == 1
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
elseif numColors == 2
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
end

% Get the ROIs
if numColors == 1
    num_glomeruli = 18;
    A = mean(stackMaxInt,3);
    Rim = zeros(size(A,1),size(A,2),3);
    Rim(:,:,2) = 2*(A-min(min(A)))./(max(max(A))-min(min(A)));
    
    % Specify the ROIs
    hf = figure;
    pgons = {};
    imshow(Rim);
    hold;
    for pgs = 1:num_glomeruli
        pgons{pgs} = impoly(gca,[10*pgs,100; 10*pgs+10, 100; 10*pgs+10, 120; 10*pgs, 120]);
    end
    wait(pgons{pgs});

    pgonPos = {};
    for pgs = 1:num_glomeruli
        pgonPos{pgs} = getPosition(pgons{pgs});
    end
    delete(hf);
elseif numColors == 2
    num_glomeruli = 18;
    A = mean(RstackMean,3);
    Rim = zeros(size(A,1),size(A,2),3);
    Rim(:,:,1) = 2*(A-min(min(A)))./(max(max(A))-min(min(A)));
    B = mean(GstackMean,3);
    Gim = zeros(size(B,1),size(B,2),3);
    Gim(:,:,2) = (B-min(min(B)))./(max(max(B))-min(min(B)));
    overlayIm = Rim+Gim;

    % Specify the ROIs
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
        pgons{pgs} = impoly(gca,[10*pgs,100; 10*pgs+10, 100; 10*pgs+10, 120; 10*pgs, 120]);
    end
    wait(pgons{pgs});

    pgonPos = {};
    for pgs = 1:num_glomeruli
        pgonPos{pgs} = getPosition(pgons{pgs});
    end
    delete(hf);
end

save(strcat(imageFilename(1:5),'ReferenceROIs'),'pgonPos');

end

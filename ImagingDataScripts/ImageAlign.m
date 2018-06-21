%% Load middle plane
% Get the 2P info
[imageFilename,imagePathname] = uigetfile('*.tif','Select the Two Photon Data');
fullpath = strcat(imagePathname,imageFilename);
info = imfinfo(fullpath);
cd(imagePathname);

numFiles = 1;
fullpath = {};
num_images = 0;
num_images_ind = zeros(numFiles,1);
info = {};
for fID = 1:numFiles
    fullpath{fID} = strcat(imagePathname,imageFilename(1:end-5),num2str(fID),'.tif');
    info{fID} = imfinfo(fullpath{fID});
    num_images_ind(fID) = numel(info{fID});
    num_images = num_images + num_images_ind(fID);
end
recSpecs = info{1}(1).ImageDescription;
planeLoc = strfind(recSpecs, 'numSlices');
num_planes = 6;

% Read in the stack
% Get the file info and # of planes and frames
width = info{1}(1).Width;
height = info{1}(1).Height;
numFrames = floor(num_images/num_planes);
num_images_ind(fID) = num_images_ind(fID) - (num_images - num_planes*numFrames); 

% Find the mean volume projection
stackMean = double(zeros(height,width,200));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(incIm+imOffset,6)==2
            stackMean(:,:,(incIm+imOffset+4)/num_planes) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

%% Perform image registration
% Create a blank stack for the registered images
stackReg = zeros(size(stackMean));

% Specify the registration image and the area to be registered
ref_im = stackMean(:,:,1);
orig_res = [0 0 1024 1024];
square_EB = [75 75 850 850];

% relative offset of position of subimages
rect_offset = [(square_EB(1)-orig_res(1))
               (square_EB(2)-orig_res(2))];

h = waitbar(0.0,'Registering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');

% Register each image
for imNow = 1:size(stackMean,3)
    
    waitbar(imNow/size(stackMean,3),h,['Loading frame# ' num2str(imNow) ' out of ' num2str(size(stackMean,3))]);
        
    % Choose subregion of image
    current_im = stackMean(:,:,imNow);
    sub_EB = imcrop(current_im,square_EB);

    % Do normalized cross correlation
    c = normxcorr2(sub_EB,ref_im);

    % offset found by correlation
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(sub_EB,2))
                   (ypeak-size(sub_EB,1))];

    % total offset
    offset = corr_offset + rect_offset;
    xoffset = offset(1);
    yoffset = offset(2);

    % coordinates for placing it in the registered stack
    xbegin = round(xoffset+1);
    xend   = round(xoffset+ size(sub_EB,2));
    ybegin = round(yoffset+1);
    yend   = round(yoffset+size(sub_EB,1));

    % put it in the stack
    stackReg(ybegin:yend,xbegin:xend,imNow) = sub_EB;
end

%% Look at the mean projection and the average projection
figure;
imshow(mean(stackReg(120:end,120:end,:),3),[min(min(mean(stackReg(75:end,75:end,:),3))) max(max(mean(stackReg(75:end,75:end,:),3)))]);
imwrite(mean(stackReg(120:end,120:end,:),3)./max(max(mean(stackReg(75:end,75:end,:),3))),'ForJose.jpg');
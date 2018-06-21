%% Did you remember to set the NO frames?

%% Load the stacks
close all;
clear
clc;

% Get the 2P info
[imageFilename,imagePathname] = uigetfile('*.tif','Select the Two Photon Data');
fullpath = strcat(imagePathname,imageFilename);
info = imfinfo(fullpath);
cd(imagePathname);

numFiles = 2; %input('Number of tifs?');
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
planeLoc = strfind(recSpecs, 'numFramesPerVolume');
discardLoc = strfind(recSpecs, 'numDiscardFlybackFrames');
num_planes = str2num(recSpecs(planeLoc+21:planeLoc+22));
num_discards = str2num(recSpecs(discardLoc+26));

% Read in the stack
% Get the file info and # of planes and frames
width = info{1}(1).Width;
height = info{1}(1).Height;
numFrames = floor(num_images/num_planes);
num_images_ind(end) = num_images_ind(end) - (num_images - num_planes*numFrames); 


% Find the mean volume projection
RstackMean = double(zeros(height,width,numFrames/2));
GstackMean = double(zeros(height,width,numFrames/2));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(ceil((incIm+imOffset)/2),num_planes) ~= 0
            if mod(incIm+imOffset,2)
                RstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}))+ RstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes)));
            else
                GstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}))+ GstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes)));
            end
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

% Find the max intensity projection
RstackMaxInt = double(zeros(height,width,numFrames/2));
GstackMaxInt = double(zeros(height,width,numFrames/2));
RminiStack = double(zeros(height,width,(num_planes-num_discards)));
GminiStack = double(zeros(height,width,(num_planes-num_discards)));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(ceil((incIm+imOffset)/2),num_planes) ~= 0
            if mod(incIm+imOffset,2)
                RminiStack(:,:,mod(ceil((incIm+imOffset)/2),num_planes)) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            else
                GminiStack(:,:,mod(ceil((incIm+imOffset)/2),num_planes)) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            end
            if (mod((incIm+imOffset)/2,num_planes) == num_planes-num_discards)
                RstackMaxInt(:,:,ceil((incIm+imOffset)/(2*num_planes))) = max(RminiStack,[],3);
                GstackMaxInt(:,:,ceil((incIm+imOffset)/(2*num_planes))) = max(GminiStack,[],3);
            end
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

%% Correct for movement, Gaussian filter, and get the raw traces for the ROIs

% Correct for movement
RstackRegMaxInt = imRegSimple(RstackMaxInt, 7);
GstackRegMaxInt = imRegSimple(GstackMaxInt, 7);
RstackRegMean = imRegSimple(RstackMean, 7);
GstackRegMean = imRegSimple(GstackMean, 7);

% Gaussian Filter
gaussianSize = [2 2];
gaussianSigma = 0.5;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
GstackXYfiltMax = double(zeros(size(GstackMaxInt)));

RstackXYfiltMean = double(zeros(size(RstackMean)));
GstackXYfiltMean = double(zeros(size(GstackMean)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackMaxInt,3)
    if mod(i,100)==0
        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
    end
    RstackXYfiltMax(:,:,i) = imfilter(RstackRegMaxInt(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMax(:,:,i) = imfilter(GstackRegMaxInt(:,:,i),Gxy,'replicate','conv');
    RstackXYfiltMean(:,:,i) = imfilter(RstackRegMean(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMean(:,:,i) = imfilter(GstackRegMean(:,:,i),Gxy,'replicate','conv');
end
delete(h);

% Background subtraction
% A = mean(RstackXYfiltMean,3);
% h = figure;
% ROIREF = roipoly(A/max(max(A)));
% delete(h);
% 
% RstackXYfiltBGsubMax = double(zeros(size(RstackXYfiltMax)));
% GstackXYfiltBGsubMax = double(zeros(size(GstackXYfiltMax)));
% RstackXYfiltBGsubMean = double(zeros(size(RstackXYfiltMean)));
% GstackXYfiltBGsubMax = double(zeros(size(GstackXYfiltMean)));
% 
% h = waitbar(0.0,'Background subtract...');
% set(h,'Position',[50 50 360 72]);
% set(h,'Name','Background subtract...');
% for i = 1:size(RstackXYfiltMax,3)
%     if mod(i,100)==0
%         waitbar(i/size(RstackXYfiltMax,3),h,['Filtering frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltMax,3))]);
%     end
%     A1 = RstackXYfiltMax(:,:,i);
%     A2 = GstackXYfiltMax(:,:,i);
%     A3 = RstackXYfiltMean(:,:,i);
%     A4 = GstackXYfiltMean(:,:,i);
%     ROI1 = A1(logical(ROIREF));
%     ROI2 = A2(logical(ROIREF));
%     ROI3 = A3(logical(ROIREF));
%     ROI4 = A4(logical(ROIREF));
%     RstackXYfiltBGsubMax(:,:,i) = RstackXYfiltMax(:,:,i) - mean2(ROI1);
%     GstackXYfiltBGsubMax(:,:,i) = GstackXYfiltMax(:,:,i) - mean2(ROI2);
%     RstackXYfiltBGsubMean(:,:,i) = RstackXYfiltMean(:,:,i) - mean2(ROI3);
%     GstackXYfiltBGsubMean(:,:,i) = GstackXYfiltMean(:,:,i) - mean2(ROI4);
% end
% delete(h);
% 
RstackXYfiltBGsubMax = RstackXYfiltMax;
GstackXYfiltBGsubMax = GstackXYfiltMax;
RstackXYfiltBGsubMean = RstackXYfiltMean;
GstackXYfiltBGsubMean = GstackXYfiltMean;

% Extract the ROIs and find the DF/F
% Fit an ellipse over the entire EB and then break it into ROI sections
A = mean(RstackXYfiltBGsubMax,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);
h = imellipse;
position = wait(h);
position(length(position)+1,:) = position(1,:);
ellipseCenter = [position(7,1) position(1,2)];
ellipseHeight = position(7,2) - position(19,2);
ellipseWidth = position(1,1) - position(13,1);
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;
delete(hf);

RROIs = zeros(num_ROIs,256,256);
for incROI = 1:num_ROIs
    xROI(1) = ellipseCenter(1);
    yROI(1) = ellipseCenter(2);
    for angs=1:5
        xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
    end
    RROIs(incROI,:,:) = roipoly(A, xROI,yROI);
end

A = mean(GstackXYfiltBGsubMax,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);
h = imellipse;
position = wait(h);
position(length(position)+1,:) = position(1,:);
ellipseCenter = [position(7,1) position(1,2)];
ellipseHeight = position(7,2) - position(19,2);
ellipseWidth = position(1,1) - position(13,1);
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;
delete(hf);

GROIs = zeros(num_ROIs,256,256);
for incROI = 1:num_ROIs
    xROI(1) = ellipseCenter(1);
    yROI(1) = ellipseCenter(2);
    for angs=1:5
        xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
    end
    GROIs(incROI,:,:) = roipoly(A, xROI,yROI);
end

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREFMax = zeros(num_ROIs,size(RstackXYfiltBGsubMax,3));
GROIaveREFMax = zeros(num_ROIs,size(GstackXYfiltBGsubMax,3));
RROIaveREFMean = zeros(num_ROIs,size(RstackXYfiltBGsubMean,3));
GROIaveREFMean = zeros(num_ROIs,size(GstackXYfiltBGsubMean,3));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(RstackXYfiltBGsubMax,3)
    if mod(i,10)==0
        waitbar(i/size(RstackXYfiltBGsubMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltBGsubMax,3))]);
    end
    for incROI = 1:num_ROIs
        AR = squeeze(RstackXYfiltBGsubMax(:,:,i));
        AG = squeeze(GstackXYfiltBGsubMax(:,:,i));
        BR = squeeze(RstackXYfiltBGsubMean(:,:,i));
        BG = squeeze(GstackXYfiltBGsubMean(:,:,i));
        RROIMax = AR(logical(squeeze(RROIs(incROI,:,:))));
        GROIMax = AG(logical(squeeze(GROIs(incROI,:,:))));
        RROIMean = BR(logical(squeeze(RROIs(incROI,:,:))));
        GROIMean = BG(logical(squeeze(GROIs(incROI,:,:))));
        RROIaveREFMax(incROI,i) = mean2(RROIMax);
        GROIaveREFMax(incROI,i) = mean2(GROIMax);
        RROIaveREFMean(incROI,i) = mean2(RROIMean);
        GROIaveREFMean(incROI,i) = mean2(GROIMean);
    end
end
delete(h);

%% Calculate deltaF over F, filter with S-G, and plot
RROIaveMax = zeros(num_ROIs,size(RstackXYfiltBGsubMax,3));
RROILowMax = RROIaveMax;
GROIaveMax = zeros(num_ROIs,size(GstackXYfiltBGsubMax,3));
GROILowMax = GROIaveMax;
RROIaveMean = zeros(num_ROIs,size(RstackXYfiltBGsubMean,3));
RROILowMean = RROIaveMean;
GROIaveMean = zeros(num_ROIs,size(GstackXYfiltBGsubMean,3));
GROILowMax = GROIaveMean;
for incROI = 1:num_ROIs
    RROILowMax = sort(squeeze(RROIaveREFMax(incROI,:)));
    GROILowMax = sort(squeeze(GROIaveREFMax(incROI,:)));
    RROILowMean = sort(squeeze(RROIaveREFMean(incROI,:)));
    GROILowMean = sort(squeeze(GROIaveREFMean(incROI,:)));
    RROIaveMax(incROI,:) = RROIaveREFMax(incROI,:)./mean(RROILowMax(1:floor(end/10)));
    GROIaveMax(incROI,:) = GROIaveREFMax(incROI,:)./mean(GROILowMax(1:floor(end/10)));
    RROIaveMean(incROI,:) = RROIaveREFMean(incROI,:)./mean(RROILowMean(1:floor(end/10)));
    GROIaveMean(incROI,:) = GROIaveREFMean(incROI,:)./mean(GROILowMean(1:floor(end/10)));
end

clear ellipse* xROI yROI

% Savitzky-Golay Filter the RFs
sgolayOrder = 3;
sgolayWindow = 7;
RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);
RROIaveMeanFilt = sgolayfilt(RROIaveMean,sgolayOrder,sgolayWindow,[],2);
GROIaveMeanFilt = sgolayfilt(GROIaveMean,sgolayOrder,sgolayWindow,[],2);

% Plot the profiles
act = figure('units','normalized','outerposition',[0 0 1 1]);

% Plot the activity
subplot(4,1,1);
imagesc(RROIaveMaxFilt);
hold on;
set(gca,'FontSize',16);
title('Red Maximum Intensity');
xlabel('Time (sec)');

subplot(4,1,2);
imagesc(RROIaveMeanFilt);
hold on;
set(gca,'FontSize',16);
title('Red Mean Intensity');
xlabel('Time (sec)');

subplot(4,1,3);
imagesc(GROIaveMaxFilt);;
hold on;
set(gca,'FontSize',16);
title('Green Maximum Intensity');
xlabel('Time (sec)');

subplot(4,1,4);
imagesc(GROIaveMeanFilt);
hold on;
set(gca,'FontSize',16);
title('Green Mean Intensity');
xlabel('Time (sec)');

colormap(brewermap(64, 'Blues'));

set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(act,strcat(imagePathname,imageFilename(1:end-4),'_BumpActNONO'),'-dpdf'); 

%% Get out the behavior information
[posFilename posPathname] = uigetfile('*.txt', 'Select the position file');
% Get the DAQ info
SYNCFilename = strcat(posFilename(1:end-4),'_SYNC',posFilename(end-3:end));
SYNCPathname = posPathname;

% Pull out the fly positional information
fileID = fopen(strcat(posPathname,posFilename));
tstamp = fgetl(fileID);
exTypeStr = strsplit(tstamp,'_');
exType = exTypeStr{end}(1:end-4);
formatSpec = '%s %f %s %f %s %f %s %f %s %f %s %f %s %d %s %d %s %d %s %d %s %d %s %f';
N=400000;
C = textscan(fileID,formatSpec,N,'CommentStyle','Current','Delimiter','\t');
t = C{1,2}; % Time
OffsetRot = C{1,4}; % Stripe rotational offset
OffsetRot = mod(OffsetRot+180, 360)-180;
OffsetFor = C{1,6}; % Stripe forward offset
OffsetLat = C{1,8}; % Stripe lateral offset
dx0 = C{1,10}; % X position of the ball from camera 1 
dx1 = C{1,12}; % X position of the ball from camera 2
dy0 = C{1,14};
dy1 = C{1,16};
closed = C{1,18};
direction = C{1,20};
trans = C{1,22};
gain = C{1,24};
fclose(fileID);

numDatPts = length(OffsetFor);

positionDat = {};
positionDat.t = t(1:numDatPts);
positionDat.OffsetRot = OffsetRot(1:numDatPts);
positionDat.OffsetFor = OffsetFor;
positionDat.OffsetLat = OffsetLat;
positionDat.dx0 = dx0;
positionDat.dx1 = dx1;
positionDat.dy0 = dy0;
positionDat.dy1 = dy1;
positionDat.closed = closed;
positionDat.direction = direction;
positionDat.trans = trans;
positionDat.gain = gain;
positionDat.exType = exType;


% Load photodiode and frame time stamps
timeStamps = importdata(strcat(SYNCPathname,SYNCFilename));

% Pull out the time stamps for the frame grab signal
clear tFrameGrab;
sampleData = 1;
upperLim = max(timeStamps(:,1));
offset = round(0.8/(88)*10000);
startFrameGrab = find(timeStamps(:,1) > 0.9*upperLim);
incDat = startFrameGrab(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,1) > 0.9*upperLim && timeStamps(incDat-2,1) < timeStamps(incDat+1,1))
        tFrameGrab(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tFrameGrab = tFrameGrab;

% Pull out the time stamps for the VR refresh
clear tVR;
sampleData = 1;
upperLim = max(timeStamps(:,2));
offset = round(0.6/(360)*10000);
VRthresh = 0.8;
startVR = find(timeStamps(:,2) > VRthresh*upperLim);
incDat = startVR(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,2) < VRthresh*upperLim && (timeStamps(incDat-1,2) < timeStamps(incDat+1,2) || timeStamps(incDat,2) < timeStamps(incDat+1,2)))
        tVR(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tVR = tVR;

% Match the behavior to the imaging
tSpan = positionDat.t-positionDat.t(1)+positionDat.tVR(1)/10000;
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFG-minFG+1,2);
OffsetForMatch = zeros(maxFG-minFG+1,1);
OffsetLatMatch = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFG+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFG+1) = positionDat.OffsetLat(tMatch(1));
end

%% Plot all of it
RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,1,1);
% imagesc(RROIaveMaxFilt);
 imagesc(RROIaveMeanFilt);
axis off;
% title('60D05');
title('37F06');
% 
subplot(4,1,2);
imagesc(GROIaveMaxFilt);
% imagesc(GROIaveMeanFilt);
axis off;
% title('37F06');
title('60D05');

subplot(4,1,3);
overlayIm = zeros([size(RROIaveMaxFilt) 3]);
% overlayIm(:,:,1) = (RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
overlayIm(:,:,1) = (RROIaveMeanFilt-min(min(RROIaveMeanFilt)))./(max(max(RROIaveMeanFilt))-min(min(RROIaveMeanFilt)));
overlayIm(:,:,2) = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
% overlayIm(:,:,2) = (GROIaveMeanFilt-min(min(GROIaveMeanFilt)))./(max(max(GROIaveMeanFilt))-min(min(GROIaveMeanFilt)));
image(overlayIm);
axis off;

colormap(brewermap(64, 'Blues'));

subplot(4,1,4);
plot(OffsetRotMatch(:,1),-OffsetRotMatch(:,2),'LineWidth',2);
axis tight;
ylabel('Global rotation (rad)','FontSize',16);
xlabel('Time (sec)','FontSize',16);
set(gca,'FontSize',16);
xlim([0 OffsetRotMatch(end,1)]);

set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(RGcomp,strcat(imagePathname,imageFilename(1:end-4),'_SummaryNONO'),'-dpdf'); 

%% For gains other than 1, convert the dx, dy values to behavioral movements

[posRot, posFor, posLat] = PositionConverter(positionDat);

% Match the position data to the framegrab times
posRotMatch = zeros(maxFG-minFG+1,1);
posForMatch = zeros(maxFG-minFG+1,1);
posLatMatch = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    posRotMatch(interp-minFG+1,2) = pi/180*posRot(tMatch(1));
    posForMatch(interp-minFG+1) = posFor(tMatch(1));
    posLatMatch(interp-minFG+1) = posLat(tMatch(1));
end

%% Plot all of it for gains other than 1
RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,1,1);
% imagesc(RROIaveMaxFilt);
 imagesc(RROIaveMeanFilt);
axis off;
% title('60D05');
title('37F06');
 
subplot(4,1,2);
imagesc(GROIaveMaxFilt);
% imagesc(GROIaveMeanFilt);
axis off;
% title('37F06');
title('60D05');

subplot(4,1,3);
overlayIm = zeros([size(RROIaveMaxFilt) 3]);
% overlayIm(:,:,1) = (RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
overlayIm(:,:,1) = (RROIaveMeanFilt-min(min(RROIaveMeanFilt)))./(max(max(RROIaveMeanFilt))-min(min(RROIaveMeanFilt)));
overlayIm(:,:,2) = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
% overlayIm(:,:,2) = (GROIaveMeanFilt-min(min(GROIaveMeanFilt)))./(max(max(GROIaveMeanFilt))-min(min(GROIaveMeanFilt)));
image(overlayIm);
axis off;

colormap(brewermap(64, 'Blues'));

subplot(4,1,4);
hold on;
plot(OffsetRotMatch(:,1),-OffsetRotMatch(:,2),'LineWidth',2,'LineStyle','--');
plot(OffsetRotMatch(:,1),-posRotMatch(:,2),'LineWidth',2);
axis tight;
ylabel('Global rotation (rad)','FontSize',16);
xlabel('Time (sec)','FontSize',16);
set(gca,'FontSize',16);
xlim([0 OffsetRotMatch(end,1)]);

set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(RGcomp,strcat(imagePathname,imageFilename(1:end-4),'_SummaryNONO'),'-dpdf'); 

%% Find PVA position and width
% Calculate and plot the circular mean
angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
angsraw = angsraw';
clear RbumpPVA;
clear GbumpPVA;
clear RbumpWidth;
clear GbumpWidth;
for ts = 1:length(RROIaveMaxFilt)
    RbumpPVA(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    GbumpPVA(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    RbumpWidth(ts) = circ_std(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    GbumpWidth(ts) = circ_std(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
end

OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
RPVAUnwrap = UnWrap(RbumpPVA,2,0);
GPVAUnwrap = UnWrap(GbumpPVA,2,0);

PVAinspect = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
hold on;
plot(OffsetRotMatch(:,1),-RPVAUnwrap(minFG:maxFG),'r','LineWidth',2);
errorbarPatch2(OffsetRotMatch(:,1),-RPVAUnwrap(minFG:maxFG), RbumpWidth(minFG:maxFG),[1 0 0]);
plot(OffsetRotMatch(:,1),-GPVAUnwrap(minFG:maxFG),'g','LineWidth',2);
errorbarPatch2(OffsetRotMatch(:,1),-GPVAUnwrap(minFG:maxFG), GbumpWidth(minFG:maxFG),[0 1 0]);
plot(OffsetRotMatch(:,1),-OffsetRotUnwrap,'k');

axis tight;
ylabel('Global rotation (rad)','FontSize',16);
xlabel('Time (sec)','FontSize',16);
set(gca,'FontSize',16);

set(PVAinspect,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PVAinspect,strcat(imagePathname,imageFilename(1:end-4),'_PVASummary'),'-dpdf'); 

%% Make a video of bump cross sections over time
writerObj = VideoWriter(strcat(imageFilename(1:end-4),'_XCNoNO.avi'));
writerObj.FrameRate= 1./mean(diff(OffsetRotMatch(:,1)));
open(writerObj);

% RXCPlot =(RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
RXCPlot =(RROIaveMeanFilt-min(min(RROIaveMeanFilt)))./(max(max(RROIaveMeanFilt))-min(min(RROIaveMeanFilt)));
GXCPlot = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
% GXCPlot = (GROIaveMeanFilt-min(min(GROIaveMeanFilt)))./(max(max(GROIaveMeanFilt))-min(min(GROIaveMeanFilt)));

XCvid = figure;

for frmNum=1:length(RROIaveMaxFilt)
    hold off;
    plot(RXCPlot(:,frmNum),'m','LineWidth',2);
    hold on;
    plot(GXCPlot(:,frmNum),'g','LineWidth',2);
    xlim([1 16]);
    ylim([0 1]);
    set(gca,'FontSize');
    axis off;
    timeNow = max(1,frmNum-minFG);
    timeNow = min(timeNow,maxFG-minFG);
    text(1, 1, strcat('t=',num2str(OffsetRotMatch(timeNow,1))),'FontSize',16);
    
    frame = getframe(XCvid);
    writeVideo(writerObj,frame);
end

delete(XCvid); 
close(writerObj);

%% Look at the different bumps for right and left turns
% Sort into right and left turn bins
RturnsG = [];
RturnsR = [];
LturnsG = [];
LturnsR = [];

velThresh = 30; % angular velocity threshold
velThresh = velThresh*pi/180*mean(diff(OffsetRotMatch(:,1)));

OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(1./mean(diff(OffsetRotMatch(:,1)))));

% RXCPlot =(RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
RXCPlot =(RROIaveMeanFilt-min(min(RROIaveMeanFilt)))./(max(max(RROIaveMeanFilt))-min(min(RROIaveMeanFilt)));
GXCPlot = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
% GXCPlot = (GROIaveMeanFilt-min(min(GROIaveMeanFilt)))./(max(max(GROIaveMeanFilt))-min(min(GROIaveMeanFilt)));

for i=minFG+1:maxFG
    if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThresh
        RturnsG = horzcat(RturnsG, GXCPlot(:,i));
        RturnsR = horzcat(RturnsR, RXCPlot(:,i));
    elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThresh
        LturnsG = horzcat(LturnsG, GXCPlot(:,i));
        LturnsR = horzcat(LturnsR, RXCPlot(:,i));
    end
end

% align the peaks
for i=1:size(RturnsG,2)
    pkPos = find(RturnsG(:,i)==max(RturnsG(:,i)));
    RturnsG(:,i) = circshift(RturnsG(:,i),num_ROIs/2-pkPos);
    RturnsR(:,i) = circshift(RturnsR(:,i),num_ROIs/2-pkPos);
end
for i=1:size(LturnsG,2)
    pkPos = find(LturnsG(:,i)==max(LturnsG(:,i)));
    LturnsG(:,i) = circshift(LturnsG(:,i),num_ROIs/2-pkPos);
    LturnsR(:,i) = circshift(LturnsR(:,i),num_ROIs/2-pkPos);
end

bumpSplits = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
hold on;
plot(mean(RturnsR,2)./max(mean(RturnsR,2)),'r')
plot(mean(RturnsG,2)./max(mean(RturnsR,2)),'g')
axis off;
title('CW');
    
subplot(2,2,3);
hold on;
plot(mean(LturnsR,2)./max(mean(LturnsR,2)),'r')
plot(mean(LturnsG,2)./max(mean(LturnsR,2)),'g')
axis off;
title('CCW');

subplot(2,2,2);
hold on;
plot(mean(LturnsR,2)./max(mean(LturnsR,2)),'--r')
plot(mean(RturnsR,2)./max(mean(RturnsR,2)),'r')
axis off;
title('37F06')
% title('60D05')

subplot(2,2,4);
hold on;
plot(mean(RturnsG,2)./max(mean(RturnsR,2)),'g')
plot(mean(LturnsG,2)./max(mean(LturnsR,2)),'--g')
axis off;
title('60D05')
% title('37F06')

set(bumpSplits,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpSplits,strcat(imagePathname,imageFilename(1:end-4),'_BumpSplits'),'-dpdf'); 

%% Look at shape of cross sections when bump amp > some threshold
RXCPlot =(RROIaveMeanFilt-min(min(RROIaveMeanFilt)))./(max(max(RROIaveMeanFilt))-min(min(RROIaveMeanFilt)));
GXCPlot = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
RXCAlign = RXCPlot;
GXCAlign = GXCPlot;

% align the peaks
for i=1:size(RXCAlign,2)
    pkPos = find(RXCAlign(:,i)==max(RXCAlign(:,i)));
    RXCAlign(:,i) = circshift(RXCAlign(:,i),round(num_ROIs/2)-pkPos);
end
for i=1:size(GXCAlign,2)
    pkPos = find(GXCAlign(:,i)==max(GXCAlign(:,i)));
    GXCAlign(:,i) = circshift(GXCAlign(:,i),round(num_ROIs/2)-pkPos);
end

bumpShape = figure('units','normalized','outerposition',[0 0 1 1]);

hold on;
plot(median(RXCAlign,2),'r')
plot(median(GXCAlign,2),'g')
axis off;
title('Median peak shapes');

set(bumpShape,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpShape,strcat(imagePathname,imageFilename(1:end-4),'_BumpShapes'),'-dpdf'); 

%% Save the data
save(strcat(imagePathname,posFilename(1:end-4),'.mat'),...
    'RROIaveMaxFilt','GROIaveMaxFilt','RROIaveMeanFilt','GROIaveMeanFilt','positionDat');

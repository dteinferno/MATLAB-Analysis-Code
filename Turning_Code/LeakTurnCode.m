%% Load the stacks
close all;
clear
clc;

% Get the 2P info
[imageFilename,imagePathname] = uigetfile('*.tif','Select the Two Photon Data');
fullpath = strcat(imagePathname,imageFilename);
info = imfinfo(fullpath);
cd(imagePathname);

numFiles = 1; %input('Number of tifs?');
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
num_discards = str2num(recSpecs(discardLoc+26))+1;

% Read in the stack
% Get the file info and # of planes and frames
width = info{1}(1).Width;
height = info{1}(1).Height;
numFrames = floor(num_images/num_planes);
num_images_ind(end) = num_images_ind(end) - (num_images - num_planes*numFrames); 


% Find the mean volume projection
Rstack = double(zeros(height,width,numFrames/2));
Gstack = double(zeros(height,width,numFrames/2));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(ceil((incIm+imOffset)/2),num_planes) == 6
            if mod(incIm+imOffset,2)
                Rstack(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            else
                Gstack(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            end
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

%% Show the two stacks and their superposition
figure;
A = squeeze(mean(Rstack,3));
A = (A-min(min(A)))./(max(max(A))-min(min(A)));
B = squeeze(mean(Gstack,3));
B = (B-min(min(B)))./(max(max(B))-min(min(B)));

subplot(1,3,1);
Rim = zeros(size(A,1), size(A,2),3);
Rim(:,:,1) = A;
imshow(Rim);

subplot(1,3,2);
Gim = zeros(size(B,1), size(B,2),3);
Gim(:,:,2) = B;
imshow(Gim);

subplot(1,3,3);
imshow(Rim+Gim)

%% Define the ROI for the two NO and calculate the intensity in the ROIs over time
num_ROIs=2

A = mean(Gstack,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);

GROIs = zeros(num_ROIs,256,256);

h1 = imellipse;
position = wait(h1);
GROIs(1,:,:) = roipoly(A, position(:,1), position(:,2));

h2 = imellipse;
position = wait(h2);
GROIs(2,:,:) = roipoly(A, position(:,1), position(:,2));
close(hf);

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREF = zeros(num_ROIs,size(Rstack,3));
GROIaveREF = zeros(num_ROIs,size(Gstack,3));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(Rstack,3)
    if mod(i,10)==0
        waitbar(i/size(Rstack,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(Rstack,3))]);
    end
    for incROI = 1:num_ROIs
        AR = squeeze(Rstack(:,:,i));
        AG = squeeze(Gstack(:,:,i));
        RROI = AR(logical(squeeze(GROIs(incROI,:,:))));
        GROI = AG(logical(squeeze(GROIs(incROI,:,:))));
        RROIaveREF(incROI,i) = mean2(RROI);
        GROIaveREF(incROI,i) = mean2(GROI);
    end
end
delete(h);


RROIave = zeros(num_ROIs,size(Rstack,3));
RROILow = RROIave;
GROIave = zeros(num_ROIs,size(Gstack,3));
GROILow = GROIave;
for incROI = 1:num_ROIs
    RROILow = sort(squeeze(RROIaveREF(incROI,:)));
    GROILow = sort(squeeze(GROIaveREF(incROI,:)));
    RROIave(incROI,:) = RROIaveREF(incROI,:)./mean(RROILow(1:floor(end/10)));
    GROIave(incROI,:) = GROIaveREF(incROI,:)./mean(GROILow(1:floor(end/10)));
end

% load position data
posFilename = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_02.txt';
% posFilename = 'Fly2_7day_6fxVT_jRGC1ax60D05_04.txt';

% Get the DAQ info
SYNCFilename = strcat(posFilename(1:end-4),'_SYNC',posFilename(end-3:end));

% Pull out the fly positional information
fileID = fopen(posFilename);
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
timeStamps = importdata(SYNCFilename);

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
num_planes = round(length(positionDat.tFrameGrab)/...
    length(RROIave));
tSpan = positionDat.t - ...
    positionDat.t(1) + ...
    positionDat.tVR(1)/10000;
minFG = ceil(min(find(positionDat.tFrameGrab >=...
    positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);

            
% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFG-minFG+1,2);
OffsetForMatch = zeros(maxFG-minFG+1,1);
OffsetLatMatch = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >=...
        (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)-...
        positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFG+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFG+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);

Rdiff = RROIave(1,minFG:maxFG)./max(RROIave(1,minFG:maxFG))-RROIave(2,minFG:maxFG)./max(RROIave(2,minFG:maxFG));
Gdiff = GROIave(1,minFG:maxFG)./max(GROIave(1,minFG:maxFG))-GROIave(2,minFG:maxFG)./max(GROIave(2,minFG:maxFG));

figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(4,1,1);
% hold on;
% plot(OffsetRotMatch(:,1),RROIave(1,minFG:maxFG)./max(RROIave(1,minFG:maxFG)),'r');
% plot(OffsetRotMatch(:,1),GROIave(1,minFG:maxFG)./max(GROIave(1,minFG:maxFG)),'g');

corrcoef(RROIave(1,minFG:maxFG),GROIave(1,minFG:maxFG))
% subplot(4,1,2);
% plot(OffsetRotMatch(1:end-1,1),diff(OffsetRotMatchUnwrap)./max(diff(OffsetRotMatchUnwrap)),'b');
% 
% subplot(4,1,3);
% hold on;
% plot(OffsetRotMatch(:,1),RROIave(2,minFG:maxFG)./max(RROIave(2,minFG:maxFG)),'r');
% plot(OffsetRotMatch(:,1),GROIave(2,minFG:maxFG)./max(GROIave(2,minFG:maxFG)),'g');

corrcoef(RROIave(2,minFG:maxFG),GROIave(2,minFG:maxFG))

subplot(2,1,1);
hold on;
plot(OffsetRotMatch(1:end-1,1),diff(OffsetRotMatchUnwrap)./max(diff(OffsetRotMatchUnwrap)),'b');
plot(OffsetRotMatch(:,1),Gdiff,'g');
set(gca,'FontSize',16);
xlim([0 65]);
ylim([-1 1]);

corrcoef(vRot,Gdiff(2:end))

subplot(2,1,2);
hold on;
plot(OffsetRotMatch(1:end-1,1),vRot./max(vRot),'b');
plot(OffsetRotMatch(:,1),Rdiff,'r');
set(gca,'FontSize',16);
xlim([0 65]);
ylim([-1 1]);
xlabel('Time (sec)');

corrcoef(vRot,Rdiff(2:end))
% subplot(4,1,3);
% hold on;
% plot(OffsetRotMatch(:,1),RROIave(1,minFG:maxFG)./max(RROIave(1,minFG:maxFG)),'r');
% plot(OffsetRotMatch(:,1),RROIave(2,minFG:maxFG)./max(RROIave(1,minFG:maxFG)),'Color',[1 0.5 0.5]);
% 
% subplot(4,1,4);
% hold on;
% plot(OffsetRotMatch(:,1),GROIave(1,minFG:maxFG)./max(GROIave(1,minFG:maxFG)),'g');
% plot(OffsetRotMatch(:,1),GROIave(2,minFG:maxFG)./max(GROIave(1,minFG:maxFG)),'Color',[0.5 1 0.5]);

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
hold on;
plot(OffsetRotMatch(1:end-1,1),smooth(diff(OffsetRotMatchUnwrap)./max(diff(OffsetRotMatchUnwrap)),round(1./mean(diff(OffsetRotMatch(:,1))))),'b');
plot(OffsetRotMatch(:,1),smooth(Gdiff,round(1./mean(diff(OffsetRotMatch(:,1))))),'g');
set(gca,'FontSize',16);
xlim([0 125]);
ylim([-1 1]);

corrcoef(smooth(vRot,round(1./mean(diff(OffsetRotMatch(:,1))))),smooth(Gdiff(2:end),round(1./mean(diff(OffsetRotMatch(:,1))))))

subplot(2,1,2);
hold on;
plot(OffsetRotMatch(1:end-1,1),smooth(vRot./max(vRot),round(1./mean(diff(OffsetRotMatch(:,1))))),'b');
plot(OffsetRotMatch(:,1),smooth(Rdiff,round(1./mean(diff(OffsetRotMatch(:,1))))),'r');
set(gca,'FontSize',16);
xlim([0 125]);
ylim([-1 1]);
xlabel('Time (sec)');

corrcoef(smooth(vRot,round(1./mean(diff(OffsetRotMatch(:,1))))),smooth(Rdiff(2:end),round(1./mean(diff(OffsetRotMatch(:,1))))))

%% Define an ROI for the red channel and calculate the intensities over time
num_ROIs = 1;

A = mean(Rstack,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);

RROImask = roipoly;

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREF2 = zeros(num_ROIs,size(Rstack,3));
GROIaveREF2 = zeros(num_ROIs,size(Gstack,3));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(Rstack,3)
    if mod(i,10)==0
        waitbar(i/size(Rstack,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(Rstack,3))]);
    end
    for incROI = 1:num_ROIs
        AR = squeeze(Rstack(:,:,i));
        AG = squeeze(Gstack(:,:,i));
        RROI = AR(logical(RROImask));
        GROI = AG(logical(RROImask));
        RROIaveREF2(incROI,i) = mean2(RROI);
        GROIaveREF2(incROI,i) = mean2(GROI);
    end
end
delete(h);


RROIave2 = zeros(num_ROIs,size(Rstack,3));
RROILow = RROIave2;
GROIave2 = zeros(num_ROIs,size(Gstack,3));
GROILow = GROIave2;
for incROI = 1:num_ROIs
    RROILow = sort(squeeze(RROIaveREF2(incROI,:)));
    GROILow = sort(squeeze(GROIaveREF2(incROI,:)));
    RROIave2(incROI,:) = RROIaveREF2(incROI,:)./mean(RROILow(1:floor(end/10)));
    GROIave2(incROI,:) = GROIaveREF2(incROI,:)./mean(GROILow(1:floor(end/10)));
end

%% Define an ROI for the green channel and calculate the intensities over time
num_ROIs = 1;

A = mean(Gstack,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);

GROImask = roipoly;

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREF3 = zeros(num_ROIs,size(Rstack,3));
GROIaveREF3 = zeros(num_ROIs,size(Gstack,3));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(Rstack,3)
    if mod(i,10)==0
        waitbar(i/size(Rstack,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(Rstack,3))]);
    end
    for incROI = 1:num_ROIs
        AR = squeeze(Rstack(:,:,i));
        AG = squeeze(Gstack(:,:,i));
        RROI = AR(logical(GROImask));
        GROI = AG(logical(GROImask));
        RROIaveREF3(incROI,i) = mean2(RROI);
        GROIaveREF3(incROI,i) = mean2(GROI);
    end
end
delete(h);


RROIave3 = zeros(num_ROIs,size(Rstack,3));
RROILow = RROIave3;
GROIave3 = zeros(num_ROIs,size(Gstack,3));
GROILow = GROIave3;
for incROI = 1:num_ROIs
    RROILow = sort(squeeze(RROIaveREF3(incROI,:)));
    GROILow = sort(squeeze(GROIaveREF3(incROI,:)));
    RROIave3(incROI,:) = RROIaveREF3(incROI,:)./mean(RROILow(1:floor(end/10)));
    GROIave3(incROI,:) = GROIaveREF3(incROI,:)./mean(GROILow(1:floor(end/10)));
end

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(OffsetRotMatch(:,1),smooth(RROIave3(1,minFG:maxFG)./max(RROIave(1,minFG:maxFG)),round(1./mean(diff(OffsetRotMatch(:,1))))),'r');
plot(OffsetRotMatch(:,1),-smooth(GROIave3(1,minFG:maxFG)./max(GROIave(1,minFG:maxFG)),round(1./mean(diff(OffsetRotMatch(:,1))))),'g');

corrcoef(smooth(RROIave3(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1))))),smooth(GROIave3(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1))))))
corrcoef(RROIave3(1,minFG:maxFG),GROIave3(1,minFG:maxFG))

%% Show the ROIs on the images
figure;
A = squeeze(mean(Rstack,3));
A = (A-min(min(A)))./(max(max(A))-min(min(A)));
B = squeeze(mean(Gstack,3));
B = (B-min(min(B)))./(max(max(B))-min(min(B)));

subplot(1,3,1);
Rim = zeros(size(A,1), size(A,2),3);
Rim(:,:,1) = A;
imshow(Rim);

subplot(1,3,2);
Gim = zeros(size(B,1), size(B,2),3);
Gim(:,:,2) = B;
imshow(Gim);

subplot(1,3,3);
imshow(Rim+Gim)
hold on;
boundpt = [0 0];
for roiTrace = 1:size(RROImask,1)
    boundpt = [0 0];
    for lineScan = 1:size(RROImask,2)
        boundptFnd = find(RROImask(:,lineScan));
        if ~isempty(boundptFnd)
            boundpt(2) = lineScan;
            boundpt(1) = boundptFnd(1);
            break;
        end
    end
    contour = bwtraceboundary(RROImask,[boundpt(1) boundpt(2)],'S');
    plot(contour(:,2),contour(:,1),'w','LineWidth',0.5);
end

for roiTrace = 1:size(GROImask,1)
    boundpt = [0 0];
    for lineScan = 1:size(GROImask,2)
        boundptFnd = find(GROImask(:,lineScan));
        if ~isempty(boundptFnd)
            boundpt(2) = lineScan;
            boundpt(1) = boundptFnd(1);
            break;
        end
    end
    contour = bwtraceboundary(GROImask,[boundpt(1) boundpt(2)],'S');
    plot(contour(:,2),contour(:,1),'w','LineWidth',0.5);
end

for roiTrace = 1:size(GROIs,1)
    boundpt = [0 0];
    for lineScan = 1:size(GROIs,2)
        boundptFnd = find(GROIs(roiTrace,:,lineScan));
        if ~isempty(boundptFnd)
            boundpt(2) = lineScan;
            boundpt(1) = boundptFnd(1);
            break;
        end
    end
    contour = bwtraceboundary(squeeze(GROIs(roiTrace,:,:)),[boundpt(1) boundpt(2)],'S');
    plot(contour(:,2),contour(:,1),'w','LineWidth',0.5);
end

%% Plot smoothed signals
Rsmooth2 = smooth(RROIaveREF2(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1)))));
Gsmooth2 = smooth(GROIaveREF2(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1)))));
Rsmooth3 = smooth(RROIaveREF3(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1)))));
Gsmooth3 = smooth(GROIaveREF3(1,minFG:maxFG),round(1./mean(diff(OffsetRotMatch(:,1)))));
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
hold on;
plot(OffsetRotMatch(:,1),(Rsmooth2-min(Rsmooth2))./(max(Rsmooth2)-min(Rsmooth2)),'r');
plot(OffsetRotMatch(:,1),(Gsmooth2-min(Gsmooth3))./(max(Gsmooth3)-min(Gsmooth3)),'g');
xlim([0 65]);
ylim([-0.4 1]);

corrcoef((Rsmooth2-min(Rsmooth2))./(max(Rsmooth2)-min(Rsmooth2)),(Gsmooth2-min(Gsmooth3))./(max(Gsmooth3)-min(Gsmooth3)))
corrcoef(RROIave2(1,minFG:maxFG),GROIave2(1,minFG:maxFG))

subplot(2,1,2)
hold on;
plot(OffsetRotMatch(:,1),(Rsmooth3-min(Rsmooth2))./(max(Rsmooth2)-min(Rsmooth2)),'r');
plot(OffsetRotMatch(:,1),(Gsmooth3-min(Gsmooth3))./(max(Gsmooth3)-min(Gsmooth3)),'g');
xlim([0 65]);
ylim([-0.4 1]);

corrcoef((Rsmooth3-min(Rsmooth2))./(max(Rsmooth2)-min(Rsmooth2)),(Gsmooth3-min(Gsmooth3))./(max(Gsmooth3)-min(Gsmooth3)))
corrcoef(RROIave3(1,minFG:maxFG),GROIave3(1,minFG:maxFG))
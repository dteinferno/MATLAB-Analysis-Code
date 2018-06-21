%% Load the stacks and the ROI data
clear;
clc;

% tifName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_00002_00001.tif';
% moveName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_02.TXT';
% allPathname = 'D:\Imaging\2Color\20160429\';
% fileName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_00002.mat';
% aviName = 'Fly2_Dark_02.avi';

tifName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003_00001.tif';
moveName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_03.TXT';
allPathname = 'D:\Imaging\2Color\20160505\';
fileName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003.mat';
aviName = 'Fly2_Dark_03.avi';

tifChunks = 1;
[fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = ...
    VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,1);
load(strcat(allPathname,fileName));
load(strcat(allPathname,fileName(1:end-4),'_ROIs.mat'));

%% Correct for movement
RstackRegMaxInt = imRegSimple(RstackMaxInt, 7);
GstackRegMaxInt = imRegSimple(GstackMaxInt, 7);
RstackRegMean = imRegSimple(RstackMean, 7);
GstackRegMean = imRegSimple(GstackMean, 7);

%% Gaussian Filter

% Gaussian Filter
gaussianSize = [5 5];
gaussianSigma = 2;
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
    RstackXYfiltMax(:,:,i) = imfilter(RstackMaxInt(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMax(:,:,i) = imfilter(GstackMaxInt(:,:,i),Gxy,'replicate','conv');
    RstackXYfiltMean(:,:,i) = imfilter(RstackMean(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMean(:,:,i) = imfilter(GstackMean(:,:,i),Gxy,'replicate','conv');
end
delete(h);

%% Make the movie
tStart = 5;
tStop = 35;
num_planes = length(positionDat.tFrameGrab)./size(GstackXYfiltMax,3);

% Open the movie object and make the figure
writerObj = VideoWriter(strcat(moveName(1:end-4),'Fixed.avi'));
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
writerObj.FrameRate= framerate;
open(writerObj);
movieFig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','k');

% Find the regions where both the activity and behavior are recorded
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
minFG = max(minFG,round(tStart*framerate));
maxFG = round(length(positionDat.tFrameGrab)/num_planes);
maxFG = min(maxFG,round(tStop*framerate));

% Create the colored images from the two stacks
RmaxVal = max(max(max(RstackXYfiltMax)));
RminVal = min(min(min(RstackXYfiltMax)));
GmaxVal = max(max(max(GstackXYfiltMax)));
GminVal = min(min(min(GstackXYfiltMax)));

curRFrame = (RstackXYfiltMax(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
curGFrame = 6*(GstackXYfiltMax(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;
overlayIm = zeros([size(curRFrame) 3]);
overlayIm(:,:,1) = curRFrame;
overlayIm(:,:,2) = curGFrame;

% Create the base figure;
subplot(4,3,1);
imshow(RIm);
axis off;
text(0,size(RIm,1)-20,'60D05','Color','w','FontSize',14);
title('Red fluorescence','Color','w','FontSize',14);
hold on;

subplot(4,3,2);
imshow(GIm);
axis off;
text(0,size(GIm,1)-20,'37F06','Color','w','FontSize',14);
title('Green fluorescence','Color','w','FontSize',14);
hold on;

subplot(4,3,3);
imshow(overlayIm);
axis off;
title('Combined fluorescence','Color','w','FontSize',14);

% Look at the cross sections
subplot(4,3,4:9);
hold on;
plot(linspace(-pi,pi,16),GROIaveMax(:,minFG)./max(max(GROIaveMax)),'g','LineWidth',2);
plot(linspace(-pi,pi,16),RROIaveMax(:,minFG)./max(max(RROIaveMax)),'r','LineWidth',2);
ylim([0 1]);
xlim([-pi pi]);
axis off;
title('\DeltaF/F for the ROIs','Color','w','FontSize',14);

% Look at the fly on the ball
subplot(4,3,11);
walkVid = VideoReader(aviName);
vidWidth = walkVid.Width;
vidHeight = walkVid.Height;
mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
mov(1).cdata = readFrame(walkVid);
for advmovie = 2:floor(minFG*num_planes/2)-1
    mov(advmovie).cdata = readFrame(walkVid);
end
imshow(mov(advmovie).cdata(50:300,100:385,1));

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFG-minFG+1,2);
OffsetForMatch = zeros(maxFG-minFG+1,1);
OffsetLatMatch = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFG+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFG+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);

% Plot the rotational velocity
subplot(4,3,10);
if vRot(1) < 0
    bar(-vRot(1),'w');
else
    bar(0,'w');
end
ylim([0 max(abs(vRot))]);
axis off;
title('v_{cw}','Color','w');

subplot(4,3,12);
if vRot(1) > 0
    bar(vRot(1),'w');
else
    bar(0,'w');
end
ylim([0 max(abs(vRot))]);
axis off;
title('v_{ccw}','Color','w');



h = waitbar(0.0,'Making movie...');
set(h,'Position',[800 50 360 72]);
set(h,'Name','Making movie...');

frame = getframe(movieFig);
writeVideo(writerObj,frame);

% Hold steady
tHold = 3;
for frameHold = 1:round(tHold*framerate)
    writeVideo(writerObj,frame);
end

% Add the ROIs
tAdd = 12;
for frameAdd = 1:round(tAdd*framerate)
    subplot(4,3,1);
    cla;
    imshow(RIm);
    axis off;
    text(0,size(RIm,1)-20,'60D05','Color','w','FontSize',14);
    title('Red fluorescence','Color','w','FontSize',14);
    hold on;
    boundpt = [0 0];
    for roiTrace = 1:size(RROIs,1)
        boundpt = [0 0];
        for lineScan = 1:size(RROIs,2)
            boundptFnd = find(RROIs(roiTrace,:,lineScan));
            if ~isempty(boundptFnd)
                boundpt(2) = lineScan;
                boundpt(1) = boundptFnd(1);
                break;
            end
        end
        contour = bwtraceboundary(squeeze(RROIs(roiTrace,:,:)),[boundpt(1) boundpt(2)],'S');
        plot(contour(:,2),contour(:,1),'w','LineWidth',0.5);
        if roiTrace == ceil(size(RROIs,1)*frameAdd/round(tAdd*framerate))
            fill(contour(:,2),contour(:,1),'w');
        end 
    end
    
    subplot(4,3,2);
    cla;
    imshow(GIm);
    axis off;
    text(0,size(GIm,1)-20,'37F06','Color','w','FontSize',14);
    title('Green fluorescence','Color','w','FontSize',14);
    hold on;
    boundpt = [0 0];
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
        if roiTrace == ceil(size(GROIs,1)*frameAdd/round(tAdd*framerate))
            fill(contour(:,2),contour(:,1),'w');
        end 
    end
    frame = getframe(movieFig);
    writeVideo(writerObj,frame);
    
    subplot(4,3,4:9)
    cla;
    hold on;
    plot(linspace(-pi,pi,16),GROIaveMax(:,minFG)./max(max(GROIaveMax)),'g','LineWidth',2);
    plot(linspace(-pi,pi,16),RROIaveMax(:,minFG)./max(max(RROIaveMax)),'r','LineWidth',2);
    ylim([0 1]);
    xlim([-pi pi]);
    axis off;
    title('\DeltaF/F for the ROIs','Color','w','FontSize',14);
    line([-pi+(ceil(size(GROIs,1)*frameAdd/round(tAdd*framerate))-1)*2*pi/15 -pi+(ceil(size(GROIs,1)*frameAdd/round(tAdd*framerate))-1)*2*pi/15], ...
        [0 max(max(max(RROIaveMax)),max(max(GROIaveMax)))],'Color','w','LineStyle','--');
end

for frameNum = minFG+1:maxFG-1
    
    waitbar((frameNum-minFG)/(maxFG-minFG-1),h,['Loading frame# ' num2str(frameNum-minFG) ' out of ' num2str(maxFG-minFG-1)]);
    
    curRFrame = (RstackXYfiltMax(:,:,frameNum)-RminVal)./(RmaxVal-RminVal);
    curGFrame = 6*(GstackXYfiltMax(:,:,frameNum)-GminVal)./(GmaxVal-GminVal);

    RIm = zeros([size(curRFrame) 3]);
    RIm(:,:,1) = curRFrame;
    GIm = zeros([size(curGFrame) 3]);
    GIm(:,:,2) = curGFrame;
    overlayIm = zeros([size(curRFrame) 3]);
    overlayIm(:,:,1) = curRFrame;
    overlayIm(:,:,2) = curGFrame;

   % Create the base figure;
    subplot(4,3,1);
    hold off;
    cla;
    imshow(RIm);
    axis off;

    subplot(4,3,2);
    hold off;
    cla;
    imshow(GIm);
    axis off;

    subplot(4,3,3);
    cla;
    imshow(overlayIm);
    axis off;

    % Look at the cross sections
    subplot(4,3,4:9);
    title '';
    cla;
    hold on;
    plot(linspace(-pi,pi,16),GROIaveMax(:,frameNum),'g','LineWidth',2);
    plot(linspace(-pi,pi,16),RROIaveMax(:,frameNum),'r','LineWidth',2);
    ylim([0 max(max(max(RROIaveMax)),max(max(GROIaveMax)))]);
    xlim([-pi pi]);
    axis off;

    % Look at the fly on the ball
    subplot(4,3,11);
    mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);
    for advmovie = floor((frameNum-1)*num_planes/2):floor(frameNum*num_planes/2)-1
        mov(advmovie).cdata = readFrame(walkVid);
    end
    imshow(mov(advmovie).cdata(50:300,100:385,1));

    % Plot the rotational velocity
    subplot(4,3,10);
    if vRot(frameNum-minFG+1) < 0
        bar(-vRot(frameNum-minFG+1),'w');
    else
        bar(0,'w');
    end
    ylim([0 max(abs(vRot))]);
    axis off;

    subplot(4,3,12);
    if vRot(frameNum-minFG+1) > 0
        bar(vRot(frameNum-minFG+1),'w');
    else
        bar(0,'w');
    end
    ylim([0 max(abs(vRot))]);
    axis off;

    frame = getframe(movieFig);
    writeVideo(writerObj,frame);
end
delete(h);
delete(movieFig);
    
close(writerObj);

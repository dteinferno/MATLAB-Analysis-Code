% 08/5/2015
% Code to look at 38C04 neurons
% Dan Turner-Evans

%% Clear out old data
clear;
clc;
close all;
    
%% Specify filenames
cd('D:/Imaging/38C04');
[posFilename posPathname] = uigetfile('*.txt', 'Select the position file');
[SYNCFilename,SYNCPathname] = uigetfile('*.txt', 'Select the SYNC file');
[imageFilename,imagePathname] = uigetfile('*.tif','Select the Two Photon Data');
numFiles = input('Number of tifs?');
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
num_planes = str2num(recSpecs(planeLoc+12));

% Pull out and plot the fly positional information
% Load the position information
fileID = fopen(strcat(posPathname,posFilename));
tstamp = fgetl(fileID)
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

%% Analyze the data!
openLoop = 0;
% Pull out timestamps for display and frame grab

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

% Pull out the time stamps for the VR refresh
clear tVR;
sampleData = 1;
upperLim = max(timeStamps(:,2));
offset = round(0.6/(360)*10000);
startVR = find(timeStamps(:,2) > 0.85*upperLim);
incDat = startVR(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,2) < 0.85*upperLim && (timeStamps(incDat-1,2) < timeStamps(incDat+1,2) || timeStamps(incDat,2) < timeStamps(incDat+1,2)))
        tVR(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

figure;
plot(timeStamps(:,2));
xlim([tVR(1)-200 tVR(1)+800]);
hold on;
scatter(tVR,upperLim*0.85+zeros(length(tVR),1),'r');

% Discriminate between ther different trials
trialbreaks1 = find(abs(diff(closed))>=1);
trialbreaks2 = find(abs(diff(direction))>=1);
trialbreaks = union(trialbreaks1,trialbreaks2);
trialbreaks = vertcat(0,trialbreaks,length(t));
numTrials = length(trialbreaks)-1;

flyPos = figure;set(flyPos,'Position',[50 50 750 900]);

% Plot the fly's position in the arena throughout the trial
subplot(5,1,1:3);
hold on;
trialSpan = trialbreaks(2)+1:trialbreaks(3);
c = zeros(length(OffsetLat(1:10:end)),1);
c(round(trialbreaks(2)/10):round(trialbreaks(3)/10)) = linspace(1, 10, round(trialbreaks(3)/10)-round(trialbreaks(2)/10)+1);
scatter(-OffsetLat(1:10:end),OffsetFor(1:10:end),1,c);
axis equal;
%rectangle('Position',[10,-2,1,4],'FaceColor','b');
viscircles([0 0], 1,'EdgeColor','b');
title('Position Throughout the Trial');

% Plot the fly's angle throughout the trial
subplot(5,1,4);
scatter(t(1:10:end),pi/180*OffsetRot(1:10:end),1,c);
xlabel('Time (sec)');
ylabel('Angular Orientation (rad)');
ylim([-pi pi]);
axis tight;
line([t(trialbreaks(2)) t(trialbreaks(2))],[-pi pi],'Color',[0 0 0]);
line([t(trialbreaks(3)) t(trialbreaks(3))],[-pi pi],'Color',[0 0 0]);

% Plot a histogram of the fly's average speed
frameRate = 120;
netSpeed = sqrt(diff(OffsetFor).^2+diff(OffsetLat).^2)*frameRate;
subplot(5,1,5);
hist(smooth(netSpeed,120),120);
title('Net Speed (per 1 sec)');
xlabel('cm/sec');
ylabel('Counts');
axis tight;


set(flyPos,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(flyPos,strcat(posPathname,posFilename(1:end-4),'_Position'),'-dpdf'); 


% Read in the stack
% Get the file info and # of planes and frames
width = info{1}(1).Width;
height = info{1}(1).Height;
numFrames = num_images/num_planes;

% Load the data
stack = double(zeros(height,width,numFrames));
h = waitbar(0.0,'Loading TIFF stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        stack(:,:,ceil((incIm+imOffset)/num_planes)) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}))+ stack(:,:,ceil((incIm+imOffset)/num_planes));
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

% Gaussian Filter, Background Subtraction, and Savitzky-Golay Filter
% Gaussian Filter
gaussianSize = [2 2];
gaussianSigma = 0.5;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
stackXYfilt = double(zeros(height,width,numFrames));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:numFrames
    if mod(i,100)==0
        waitbar(i/numFrames,h,['Filtering frame# ' num2str(i) ' out of ' num2str(numFrames)]);
    end
    stackXYfilt(:,:,i) = imfilter(stack(:,:,i),Gxy,'replicate','conv');
end
delete(h);

% Background subtraction
A = mean(stackXYfilt,3);
h = figure;
ROIREF = roipoly(A/max(max(A)));
delete(h);

stackXYfiltBGsub = double(zeros(height,width,numFrames));

h = waitbar(0.0,'Background subtract...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Background subtract...');
for i = 1:numFrames
    if mod(i,100)==0
        waitbar(i/numFrames,h,['Filtering frame# ' num2str(i) ' out of ' num2str(numFrames)]);
    end
    A = stackXYfilt(:,:,i);
    ROI = A(logical(ROIREF));
    stackXYfiltBGsub(:,:,i) = stackXYfilt(:,:,i) - mean2(ROI);
end
delete(h);

% Savitzky-Golay Filter
sgolayOrder = 3;
sgolayWindow = 7;
stackXYTfiltBGsub = sgolayfilt(stackXYfiltBGsub,sgolayOrder,sgolayWindow,[],3);

%% Extract the ROIs and find the DF/F
% Show the mean intensity and place some ROIs
A = mean(stackXYTfiltBGsub,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);

num_ellipses = input('Number of ellipses?');
num_pgons = input('Number of polygons?')
num_ROIs = num_ellipses+num_pgons;

ROIs = zeros(num_ROIs,256,256);
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
    
% Plot the ROI on each figure and calculate the average
ROIaveREF = zeros(num_ROIs,numFrames);
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:numFrames
    if mod(i,100)==0
        waitbar(i/numFrames,h,['Calculating frame# ' num2str(i) ' out of ' num2str(numFrames)]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(stackXYTfiltBGsub(:,:,i));
        ROI = A(logical(squeeze(ROIs(incROI,:,:))));
        ROIaveREF(incROI,i) = mean2(ROI);
    end
end
delete(h);
ROIave = zeros(num_ROIs,numFrames);
for incROI = 1:num_ROIs
    ROILow = sort(squeeze(ROIaveREF(incROI,:)));
    ROIave(incROI,:) = ROIaveREF(incROI,:)./sum(ROILow(1:floor(end/10)));
end

% Plot the intensity over time
figure('units','normalized','outerposition',[0 0 1 1])
for roiplots = 1:num_ROIs
    subplot(num_ROIs+1,4,[1:3]+4*(roiplots-1));
    plot(tFrameGrab(1:num_planes:end)/10000,ROIave(roiplots,:));
    axis tight;
    subplot(num_ROIs+1,4,4*roiplots);
    contourf(flipud(squeeze(ROIs(roiplots,:,:))));
    axis equal;
    axis off;
end

%% For open loop stimuli, find the fly's behavior
openLoop = 1;
Cam1RotCalibfact = 0.88;
Cam2RotCalibfact = 0.94;
Cam1PosCalibfact = 162;
Cam2PosCalibfact = 152;

flyAng = 30 * pi / 180;
dxmod0 = double(dx0) .* cos(flyAng) + double(dy0) .* sin(flyAng);
dymod0 = double(dy0) .* cos(flyAng) + double(dx0) .* sin(flyAng);
dxmod1 = double(dx1) .* cos(flyAng) + double(dy1).* sin(-flyAng);
dymod1 = double(dy1) .* cos(flyAng) + double(dx1) .* sin(-flyAng);
deltaFor = (dymod0 / Cam1PosCalibfact + dymod1 / Cam2PosCalibfact)*sqrt(2) / 2;
deltaSide = (dymod0 / Cam1PosCalibfact - dymod1 / Cam2PosCalibfact)*sqrt(2) / 2;
BallOffsetRot(1) = OffsetRot(1);
RotOffset = 2*OffsetRot(1) - (dxmod0(1) * Cam1RotCalibfact - dxmod1(1) * Cam2RotCalibfact) / 2;
BallOffsetFor(1) = deltaFor(1)*cos(BallOffsetRot(1) * pi / 180) + deltaSide(1)*sin(BallOffsetRot(1) * pi / 180);
BallOffsetSide(1) = deltaFor(1)*sin(BallOffsetRot(1) * pi / 180) - deltaSide(1)*cos(BallOffsetRot(1) * pi / 180);
for i = 2:length(dxmod0)
    BallOffsetRot(i) = (dxmod0(i) * Cam1RotCalibfact + dxmod1(i) * Cam2RotCalibfact) / 2 + RotOffset; 
    BallOffsetFor(i) = (deltaFor(i)-deltaFor(i-1))*cos(BallOffsetRot(i) * pi / 180) + (deltaSide(i)-deltaSide(i-1))*sin(BallOffsetRot(i) * pi / 180)+BallOffsetFor(i-1);
    BallOffsetSide(i) = (deltaFor(i)-deltaFor(i-1))*sin(BallOffsetRot(i) * pi / 180) - (deltaSide(i)-deltaSide(i-1))*cos(BallOffsetRot(i) * pi / 180)+BallOffsetSide(i-1);
end
BallOffsetRot = mod(BallOffsetRot+180, 360)-180;
figure; plot(OffsetRot); hold on; plot(BallOffsetRot,'r');
figure; plot(OffsetFor); hold on; plot(-5+BallOffsetFor,'r');
figure; plot(OffsetLat); hold on; plot(BallOffsetSide,'r');

%% Plot a bunch of things for each ROI

% Smooth the global rotation
OffsetRotPlot =pi/180* smooth(OffsetRot,10);

% Find the object orientation
relAng = atan2(OffsetFor,-OffsetLat);
netAngle(find(-OffsetLat > 0)) = -relAng(find(-OffsetLat > 0))-pi/2+pi/180*OffsetRot(find(-OffsetLat > 0));
netAngle(find(-OffsetLat < 0)) = -relAng(find(-OffsetLat < 0))+3*pi/2+pi/180*OffsetRot(find(-OffsetLat < 0));
netAngle=mod(netAngle+pi,2*pi)-pi;
netAngle = smooth(netAngle,10);

% Get rid of wrap-around artifacts
for i=1:length(netAngle)-1
     if abs(netAngle(i)-netAngle(i+1)) > pi/8
        netAngle(i) = NaN;
     end
     if abs(OffsetRotPlot(i)-OffsetRotPlot(i+1)) > pi
        OffsetRotPlot(i) = NaN;
     end
end

tSpan = t-t(1)+(tVR(1)-tFrameGrab(1))/10000;

% Calculate the rotational and forward velocities
if (openLoop)
   DeltaRot = diff(BallOffsetRot);
else
    DeltaRot = diff(OffsetRot);
end
notHuge = find(abs(DeltaRot) < 5);
OffsetRotDPlot = abs(smooth(DeltaRot(notHuge)*frameRate,360));

OffsetForSmooth = smooth(OffsetFor,10);
OffsetLatSmooth = smooth(OffsetLat,10);
netSpeed = sqrt(diff(OffsetForSmooth).^2+diff(OffsetLatSmooth).^2)*frameRate;

for ROINum = 1:num_ROIs
    ROINow = figure('units','normalized','outerposition',[0 0 1 1])
    
    % Plot the fly's trajectory
    subplot(4,3,1);
    hold on;
    c = zeros(length(OffsetLat(1:10:end)),1);
    c(round(trialbreaks(2)/10):round(trialbreaks(3)/10)) = linspace(1, 10, round(trialbreaks(3)/10)-round(trialbreaks(2)/10)+1);
    scatter(-OffsetLat(1:10:end),OffsetFor(1:10:end),1,c);
    axis equal;
    viscircles([0 0], 1,'EdgeColor','b');
    title('Position Throughout the Trial');

    % Plot the image
    subplot(4,3,4);
    A = mean(stackXYTfiltBGsub,3);
    imshow(A,[0 max(max(A))]);
    
    % Plot the ROI
    subplot(4,3,7);
    contourf(flipud(squeeze(ROIs(ROINum,:,:))));
    colormap('parula')
    axis equal;
    axis off;
    
    % Plot the ROI intensity throughout the trial and the color coded
    % position
    subplot(4,3,2:3);
    plot(tFrameGrab(1:num_planes:end)/10000,ROIave(roiplots,:),'g');
    axis tight;
    set(gca,'FontSize',16);
    hold on;
    tprog = zeros(1,length(c))+min(ROIave(roiplots,:));
    scatter(tSpan(1:10:end),tprog,3,c);
    
    % Plot the fly's global and local angle throughout the trial
    subplot(4,3,5:6);
    hold on;
    if (openLoop)
        plotyy(tSpan,pi/180*BallOffsetRot,tSpan,BallOffsetFor);
        xlim([tSpan(1),tSpan(end)]);
        ylim([-pi pi]);
        set(gca,'FontSize',16);
        legend('Global Angle','Forward Motion');
    else
        plot(tSpan,OffsetRotPlot,'g');
        plot(tSpan(1:length(netAngle)),netAngle,'r');
        xlim([tSpan(1),tSpan(end)]);
        ylim([-pi pi]);
        set(gca,'FontSize',16);
        legend('Global Angle','Object Angle');
    end
    
    % Plot the fly's rotational velocity and translational speed
    subplot(4,3,8:9);
    hold on;
    plot(tSpan(1:length(OffsetRotDPlot)),OffsetRotDPlot./max(OffsetRotDPlot),'k');
    plot(tSpan(1:length(netSpeed)),smooth(netSpeed,360)./max(smooth(netSpeed,360)),'m');
    xlim([tSpan(1),tSpan(end)]);
    ylim([0 1]);
    set(gca,'FontSize',16);
    legend('Normalized Rotational Velocity','Normalized Speed');
    xlabel('Time');
    
    % Plot the open loop values
    subplot(4,3,11:12);
    if (openLoop)
       plotyy(tSpan,t.*360./gain.*(1-double(trans)).*double(direction),tSpan,t.*gain.*double(trans).*double(direction))
       xlim([tSpan(1),tSpan(end)]);
       set(gca,'FontSize',16);
       legend('Rotation Angle','Forward Translation');
    end
    
    set(ROINow,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(ROINow,strcat(posPathname,posFilename(1:end-4),'_Activity_',num2str(ROINum)),'-dpdf'); 

    % Do correlation analysis
    t1 = tFrameGrab(1:num_planes:end)/10000;
    t2 = tSpan(1:length(OffsetRotDPlot));
    ROIforCC = interp1(t1,(ROIave(ROINum,:)-min(ROIave(ROINum,:)))/(max(ROIave(ROINum,:))-min(ROIave(ROINum,:))), t2);
    for i = 1:101
        Acoef = 0.01*(i-1);
        Bcoef = 0.01*(101-i);
        testSig = Acoef*OffsetRotDPlot./max(OffsetRotDPlot)+Bcoef*smooth(netSpeed(1:length(OffsetRotDPlot)),360)./max(smooth(netSpeed,360));
        checkVals = corrcoef(ROIforCC(~isnan(ROIforCC)),testSig(~isnan(ROIforCC)));
        findMaxCC(i) = checkVals(2,1);
    end
    bestVal = find(findMaxCC == max(findMaxCC));
    Acoef = 0.01*(bestVal-1);
    Bcoef = 0.01*(101-bestVal);
    AmpComp = figure;
    plot(tFrameGrab(1:num_planes:end)/10000,(ROIave(ROINum,:)-min(ROIave(ROINum,:)))/(max(ROIave(ROINum,:))-min(ROIave(ROINum,:))),'g');
    axis tight;
    set(gca,'FontSize',16);
    hold on;
    plot(tSpan(1:length(OffsetRotDPlot)),Acoef*OffsetRotDPlot./max(OffsetRotDPlot)+Bcoef*smooth(netSpeed(1:length(OffsetRotDPlot)),360)./max(smooth(netSpeed,360)),'b');
    legend('Activity amplitude',strcat(num2str(Acoef),'*V_{rot}+',num2str(Bcoef),'*V_{for}'));
    
    set(AmpComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(AmpComp,strcat(posPathname,posFilename(1:end-4),'_AmpComp_',num2str(ROINum)),'-dpdf'); 

end

%% Save the data
save(strcat(imagePathname,imageFilename(1:end-4),'.mat'),'num_ROIs','ROIREF','ROIs','ROIave','num_planes'); 
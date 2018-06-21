allPathname = 'D:\Imaging\2Color\20160429\';
tifName = 'Fly3_7day_6fx37F06_jRGC1ax60D05_Dark_00005_00001.tif';
moveName = 'Fly3_7day_6fx37F06_jRGC1ax60D05_Dark_05.TXT';

tifChunks = 1;

[fullpath, RstackFlat, GstackFlat, positionDat] = VRImDatLoadDiscardNoPromptOnePlane(tifName,allPathname,moveName,allPathname,tifChunks,6);
[fullpath, RstackFull, GstackFull, RstackMeanFull, GstackFullMean, positionDat] = VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,1);;
load(strcat(tifName(1:end-10),'_NOROIs','.mat'));
num_ROIs=size(ROIs,1);
num_planes = length(positionDat.tFrameGrab)/length(RstackFlat);

% Gaussian Filter
gaussianSize = [2 2];
gaussianSigma = 0.5;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RstackXYfiltFlat = double(zeros(size(RstackFlat)));
GstackXYfiltFlat = double(zeros(size(GstackFlat)));
RstackXYfiltFull = double(zeros(size(RstackFull)));
GstackXYfiltFull = double(zeros(size(GstackFull)));


h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackFlat,3)
    if mod(i,100)==0
        waitbar(i/length(RstackFlat),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackFlat))]);
    end
    RstackXYfiltFlat(:,:,i) = imfilter(RstackFlat(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltFlat(:,:,i) = imfilter(GstackFlat(:,:,i),Gxy,'replicate','conv');
    RstackXYfiltFull(:,:,i) = imfilter(RstackFull(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltFull(:,:,i) = imfilter(GstackFull(:,:,i),Gxy,'replicate','conv');
end
delete(h);

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREFFlat = zeros(num_ROIs,size(RstackXYfiltFlat,3));
GROIaveREFFlat = zeros(num_ROIs,size(GstackXYfiltFlat,3));
RROIaveREFFull = zeros(num_ROIs,size(RstackXYfiltFull,3));
GROIaveREFFull = zeros(num_ROIs,size(GstackXYfiltFull,3));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(RstackXYfiltFlat,3)
    if mod(i,10)==0
        waitbar(i/size(RstackXYfiltFlat,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltFlat,3))]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(RstackXYfiltFlat(:,:,i));
        B = squeeze(RstackXYfiltFull(:,:,i));
        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
        RROIaveREFFlat(incROI,i) = mean2(ROIMax);
        ROIMax = B(logical(squeeze(ROIs(incROI,:,:))));
        RROIaveREFFull(incROI,i) = mean2(ROIMax);

        A = squeeze(GstackXYfiltFlat(:,:,i));
        B = squeeze(GstackXYfiltFull(:,:,i));
        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
        GROIaveREFFlat(incROI,i) = mean2(ROIMax);
        ROIMax = B(logical(squeeze(ROIs(incROI,:,:))));
        GROIaveREFFull(incROI,i) = mean2(ROIMax);
    end
end
delete(h);

RROIave = zeros(num_ROIs,length(RstackFlat));
GROIave = zeros(num_ROIs,length(GstackFlat));
for incROI = 1:num_ROIs
    ROILow = sort(squeeze(RROIaveREFFull(incROI,:)));
    RROIave(incROI,:) = RROIaveREFFlat(incROI,:)./mean(abs(ROILow(1:floor(end/10))));

    ROILow = sort(squeeze(GROIaveREFFull(incROI,:)));
    GROIave(incROI,:) = GROIaveREFFlat(incROI,:)./mean(abs(ROILow(1:floor(end/10))));
end

% Do some more preprocessing
% Savitzky-Golay Filter the RFs
sgolayOrder = 3;
sgolayWindow = 11;
% RROIaveFilt = sgolayfilt(RROIave,sgolayOrder,sgolayWindow,[],2);
% GROIaveFilt = sgolayfilt(GROIave,sgolayOrder,sgolayWindow,[],2);
RROIaveFilt = sgolayfilt(RROIaveREFFlat,sgolayOrder,sgolayWindow,[],2);
GROIaveFilt = sgolayfilt(GROIaveREFFlat,sgolayOrder,sgolayWindow,[],2);

% Match the behavior to the imaging
num_planes = round(length(positionDat.tFrameGrab)/length(RROIave));
tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
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
% Unwrap the rotation
OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotUnwrap)./mean(diff(OffsetRotMatch(:,1)));
vRot = sgolayfilt(vRot,sgolayOrder,sgolayWindow);


figure;
subplot(5,1,1);
plot(OffsetRotMatch(:,1),GROIaveFilt(1,minFG:maxFG),'g');
% ylim([0 1]);
ylim([0 20]);
subplot(5,1,2);
plot(OffsetRotMatch(:,1),RROIaveFilt(1,minFG:maxFG),'r');
% ylim([0 1]);
ylim([0 80]);
subplot(5,1,3);
plot(OffsetRotMatch(:,1),GROIaveFilt(2,minFG:maxFG),'g');
% ylim([0 1]);
ylim([0 20]);
subplot(5,1,4);
plot(OffsetRotMatch(:,1),RROIaveFilt(2,minFG:maxFG),'r');
% ylim([0 1]);
ylim([0 80]);
subplot(5,1,5);
plot(OffsetRotMatch(1:end-1,1),vRot);
ylim([-3*pi/2 3*pi/2]);

%% Alternate R REF ROI
figure;
image = squeeze(mean(RstackFlat,3));
imagesc(image,[min(min(image)) max(max(image))]);
RROI = roipoly();


% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
RROIaveREFFlat2 = zeros(1,size(RstackXYfiltFlat,3));
GROIaveREFFlat2 = zeros(1,size(GstackXYfiltFlat,3));

h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:size(RstackXYfiltFlat,3)
    if mod(i,10)==0
        waitbar(i/size(RstackXYfiltFlat,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfiltFlat,3))]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(RstackXYfiltFlat(:,:,i));
        ROIMax = A(logical(squeeze(RROI)));
        RROIaveREFFlat2(i) = mean2(ROIMax);
        
        B = squeeze(GstackXYfiltFlat(:,:,i));
        ROIMax = B(logical(squeeze(RROI)));
        GROIaveREFFlat2(i) = mean2(ROIMax);
    end
end
delete(h);

RROIaveFilt2 = sgolayfilt(RROIaveREFFlat2,sgolayOrder,sgolayWindow,[],2);
GROIaveFilt2 = sgolayfilt(GROIaveREFFlat2,sgolayOrder,sgolayWindow,[],2);

figure;
subplot(5,2,1:2);
plot(OffsetRotMatch(:,1),RROIaveFilt2(minFG:maxFG),'r');
ylim([0 600]);

subplot(5,2,[3 5 7 9]);
imagesc(image,[min(min(image)) max(max(image))]);
axis equal;

subplot(5,2,[4 6 8 10]);
imagesc(RROI);
axis equal;

figure;
subplot(5,2,1:2);
plot(OffsetRotMatch(:,1),GROIaveFilt2(minFG:maxFG),'r');
ylim([0 15]);

subplot(5,2,[3 5 7 9]);
imagesc(image,[min(min(image)) max(max(image))]);
axis equal;

subplot(5,2,[4 6 8 10]);
imagesc(RROI);
axis equal;
%% Load the data

tifName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00005_00001.tif';
allPathname = 'D:\Imaging\2Color\20160907\';
moveName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_05.txt';
tifChunks = 1;
[fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,0);

fileName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00005.mat';
load(strcat(allPathname,fileName));
load(strcat(allPathname,fileName(1:4),'_ReferenceROIs.mat'));

%% Gaussian Filter
gaussianSize = [20 20];
gaussianSigma = 10;
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

%% Plot cross sections

tStart = 33.85;%11.6;
tStep =0.4;

% Find the regions where both the activity and behavior are recorded
num_planes = length(positionDat.tFrameGrab)/size(RstackMaxInt,3);
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFG = round(length(positionDat.tFrameGrab)/num_planes);

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
% vRot = tStep*framerate*diff(OffsetRotMatchUnwrap(1:round(tStep*framerate):end));
sgolayOrder = 3;
sgolayWindow = 11;
vRot = diff(OffsetRotMatchUnwrap)./mean(diff(OffsetRotMatch(:,1)));
vRot = sgolayfilt(vRot,sgolayOrder,sgolayWindow);
framerate = 1/mean(diff(OffsetRotMatch(:,1)));


% Create the colored images from the two stacks
Rmax = max(RstackXYfiltMax(:,:,minFG:maxFG),[],3);
RmaxVal = max(max(max(RstackXYfiltMax(:,:,minFG:maxFG))));
Rmin = min(RstackXYfiltMax(:,:,minFG:maxFG),[],3);
RminVal = min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
Gmax = max(GstackXYfiltMax(:,:,minFG:maxFG),[],3);
GmaxVal = max(max(max(GstackXYfiltMax(:,:,minFG:maxFG))));
Gmin = min(GstackXYfiltMax(:,:,minFG:maxFG),[],3);
GminVal = min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));

XCfig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','w');

num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

for bumpStep = 1:6
    
    frameNum = round((tStart+(bumpStep-1)*tStep)*framerate);
    curRFrame = (RstackXYfiltMax(:,:,frameNum)-Rmin)./(RmaxVal-RminVal);
    curGFrame = (GstackXYfiltMax(:,:,frameNum)-Gmin)./(GmaxVal-GminVal);

    RIm = ones([size(curRFrame) 3]);
    RIm(:,:,2) = RIm(:,:,2)-curRFrame;
    RIm(:,:,3) = RIm(:,:,3)-curRFrame;
%     RIm = zeros([size(curRFrame) 3]);
%     RIm(:,:,1) = curRFrame;
    GIm = ones([size(curGFrame) 3]);
    GIm(:,:,1) = GIm(:,:,1)-curGFrame;
    GIm(:,:,3) = GIm(:,:,3)-curGFrame;
%     GIm = zeros([size(curGFrame) 3]);
%     GIm(:,:,2) = curGFrame;

    subplot(12,4,[1+(bumpStep-1)*8 5+(bumpStep-1)*8])
%     imshow(RIm);
    imagesc(curRFrame);
    axis off;
    axis equal;
    if bumpStep == 1
        title('red fluorescence');
    end
    caxis([0 0.5]);

    subplot(12,4,[2+(bumpStep-1)*8 6+(bumpStep-1)*8])
%     imshow(GIm);
    imagesc(curGFrame);
    axis off;
    axis equal;
    if bumpStep == 1
        title('green fluorescence');
    end
    caxis([0 0.75]);
    
    subplot(12,4,[3:4]+(bumpStep-1)*8)
    hold on;
    RXC = zeros(18,1);
    RXC([2:17]) = (RROIaveMax([2:17],frameNum)-1);
    imagesc(RXC');
    line([(1+circ_mean(angsraw,RROIaveMax([2:9],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)...
        (1+circ_mean(angsraw,RROIaveMax([2:9],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    line([num_ROIs+1+(circ_mean(angsraw,RROIaveMax([10:17],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)...
        (num_ROIs+1+circ_mean(angsraw,RROIaveMax([10:17],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    
    
    caxis([0 0.75]);
    colorbar;
    axis off;
    
    subplot(12,4,[7:8]+(bumpStep-1)*8)
    hold on;
    GXC = zeros(18,1);
    GXC([1:8 11:18]) = (GROIaveMax([1:8 11:18],frameNum)-1);
    imagesc(GXC');
    line([(circ_mean(angsraw,GROIaveMax([1:8],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)...
        (circ_mean(angsraw,GROIaveMax([1:8],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    line([num_ROIs+2+(circ_mean(angsraw,GROIaveMax([11:18],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)...
        (num_ROIs+2+circ_mean(angsraw,GROIaveMax([11:18],frameNum))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    caxis([0 1.5]);
    colorbar;
%     plot(linspace(-pi,pi,16),GROIaveMax(:,frameNum)-1,'g','LineWidth',2);
%     plot(linspace(-pi,pi,16),RROIaveMax(:,frameNum)-1,'r','LineWidth',2);
%     ylim([-0.1 max(max(max(RROIaveMax)),max(max(GROIaveMax)))-1]);
%     xlim([-pi pi]);
    if bumpStep == 6
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');
    else
        axis off;
    end
    title(strcat('v_{rot} =',num2str(180/pi*vRot((frameNum-minFG)))));
end
% 
blues = colormap(brewermap(64, 'Blues'));
% greens = blues;
% greens(:,1) = blues(:,1);
% greens(:,2) = blues(:,3);
% greens(:,3) = blues(:,2);
% reds = blues;
% reds(:,1) = blues(:,3);
% reds(:,2) = blues(:,1);
% reds(:,3) = blues(:,1);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
% hot = colormap('hot');
% green_hot = hot;
% green_hot(:,1) = hot(:,2);
% green_hot(:,2) = hot(:,1);
% colormap(green_hot);
% colormap(greens);
% colormap(reds);
colormap(magentas);

% colormap('winter');
set(XCfig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(XCfig,'..\Results\PBExampleXCs_magentas','-dpdf');

%% Plot the ROIs over the image;
bumpStep = 1;
frameNum = round((tStart+(bumpStep-1)*tStep)*framerate);
curRFrame = (RstackXYfiltMax(:,:,frameNum)-Rmin)./(RmaxVal-RminVal);
curGFrame = (GstackXYfiltMax(:,:,frameNum)-Gmin)./(GmaxVal-GminVal);

test = figure;
subplot(1,2,1);
hold on;
imagesc(curRFrame);
axis off;
axis equal;
caxis([0 1]);
for pgs = 1:18
    pgons{pgs} = impoly(gca,pgonPos{pgs});
end

subplot(1,2,2);
hold on;
imagesc(curGFrame);
axis off;
axis equal;
caxis([0 1.25]);
for pgs = 1:18
    pgons{pgs} = impoly(gca,pgonPos{pgs});
end

hot = colormap('hot');
% green_hot = hot;
% green_hot(:,1) = hot(:,2);
% green_hot(:,2) = hot(:,1);
% colormap(green_hot);


set(test,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(test,'..\Results\PBExampleROIs_hot','-dpdf');
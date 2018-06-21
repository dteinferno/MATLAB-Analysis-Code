%% Load the stacks and the ROI data
% tifName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_00002_00001.tif';
% moveName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_02.TXT';
% allPathname = 'D:\Imaging\2Color\20160429\';
% fileName = 'Fly2_7day_6fx37F06_jRGC1ax60D05_Dark_00002.mat';
% aviName = 'Fly2_Dark_02.avi';

tifName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003_00001.tif';
moveName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_03.TXT';
allPathname = 'D:\Imaging\2Color\20160505\';
fileName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003.mat';
% aviName = 'Fly2_Dark_03.avi';

% tifName = 'Fly1_9day_6fx37F06_jRGC1ax60D05_Dark_00006_00001.tif';
% moveName = 'Fly1_9day_6fx37F06_jRGC1ax60D05_Dark_06.TXT';
% allPathname = 'D:\Imaging\2Color\20160423-Flip\';
% fileName = 'Fly1_9day_6fx37F06_jRGC1ax60D05_Dark_00006.mat';
% aviName = 'Fly2_Dark_02.avi';

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
circCent = [128 128];
circHeight = 100;
circWidth = 100;
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;

num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs;
angsraw = angsraw';

tStart = 47;
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
    caxis([0 0.25]);
    
    subplot(12,4,[3:4]+(bumpStep-1)*8)
    hold on;
%     imagesc((RROIaveMax(:,frameNum)-1)');
    RPlot = (RROIaveMax(:,frameNum)-1);
    for incROI = 1:num_ROIs
       xWedge(1) = circCent(1);
       yWedge(1) = circCent(2);
       for angs=1:5
           xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
       end
       patch(xWedge,yWedge,RPlot(incROI));
    end
    axis off;
    axis equal;
    caxis([0 1]);
    colorbar;
    
    subplot(12,4,[7:8]+(bumpStep-1)*8)
    hold on;
%     imagesc((GROIaveMax(:,frameNum)-1)');
    GPlot = (GROIaveMax(:,frameNum)-1);
    for incROI = 1:num_ROIs
       xWedge(1) = circCent(1);
       yWedge(1) = circCent(2);
       for angs=1:5
           xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
       end
       patch(xWedge,yWedge,GPlot(incROI));
    end
    axis off;
    axis equal;
    caxis([0 1]);
    colorbar;
    PVAang = circ_mean(angsraw, GPlot)-deltaAng/2;
    line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
    PVAang = circ_mean(angsraw, RPlot)-deltaAng/2;
    line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');
    
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
    title(strcat('v_{rot} =',num2str(180/pi*vRot((frameNum-minFG+1)))));
end

blues = colormap(brewermap(64, 'Blues'));
greens = blues;
greens(:,1) = 1-blues(:,3);
greens(:,2) = 1-blues(:,1);
greens(:,3) = 1-blues(:,3);
% reds = blues;
% reds(:,1) = blues(:,3);
% reds(:,2) = blues(:,1);
% reds(:,3) = blues(:,1);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
colormap(magentas);
% hot = colormap('hot');
% green_hot = hot;
% green_hot(:,1) = hot(:,2);
% green_hot(:,2) = hot(:,1);
% colormap(green_hot);
% colormap(greens);
% colormap(reds);
% colormap(magentas);
% colormap('winter');
set(XCfig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(XCfig,'..\Results\ExampleXCs5_magentas','-dpdf');

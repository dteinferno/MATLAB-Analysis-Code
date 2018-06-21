%% Load the zoomed out view
[ovFilename,ovPathname] = uigetfile('*.tif','Select the Two Photon Data');
fullpath = strcat(ovPathname,ovFilename);
info = imfinfo(fullpath);
cd(ovPathname);

% Load the data
width = info(1).Width;
height = info(1).Height;
numFrames = 2;

stack = double(zeros(height,width,numFrames));
for incIm = 1:numel(info)
    stack(:,:,incIm) = double(imread(fullpath, incIm, 'Info', info));
end

numEx = input('Number of excitation patterns?:');
ROIs = {};
for roi = 1:numEx
    ROIfname = uigetfile('*.mat','Select the ROI'); 
    load(ROIfname);
    ROIs{roi} = DMD_pattern_data.other_info;
end


% Create a R/G colormap
colorRat = max(max(stack(:,:,1)))/max(max(stack(:,:,2)));
map = zeros(105,3);
for git = 1:51
    map(git,1) = (git-1)/51;
    map(git+54,2) = (git-1)/51;
end
colormap(map);
axis off;

upperLim = max(max(stack(:,:,1)));
ovFig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
imshow(stack(:,:,1),[0 2*upperLim]);
% hold on;
% for roiNow=1:length(ROIs)
%     roiPts = ROIs{roiNow}.ROIs(1).drawings.position;
%     if strcmp(ROIs{roiNow}.ROIs(1).drawings.class,'impoly');
%         for roiPt = 1:length(roiPts)
%             if roiPt == length(roiPts)
%                 line([roiPts(roiPt,1) roiPts(1,1)], [roiPts(roiPt,2) roiPts(1,2)],'Color',[0 0 1]);
%             else
%                 line([roiPts(roiPt,1) roiPts(roiPt+1,1)], [roiPts(roiPt,2) roiPts(roiPt+1,2)],'Color',[0 0 1]);
%             end
%         end
%     else
%         rectangle('Position',roiPts,'Curvature',[1,1],'EdgeColor','b');
%     end
% end

subplot(1,2,2);
imshow(stack(:,:,2)*colorRat+upperLim*1.1,[0 2*upperLim])
hold on;
for roiNow=1:length(ROIs)
    roiPts = ROIs{roiNow}.ROIs(1).drawings.position;
    if strcmp(ROIs{roiNow}.ROIs(1).drawings.class,'impoly');
        for roiPt = 1:length(roiPts)
            if roiPt == length(roiPts)
                line([roiPts(roiPt,1) roiPts(1,1)], [roiPts(roiPt,2) roiPts(1,2)],'Color',[0 0 1]);
            else
                line([roiPts(roiPt,1) roiPts(roiPt+1,1)], [roiPts(roiPt,2) roiPts(roiPt+1,2)],'Color',[0 0 1]);
            end
        end
    else
        rectangle('Position',roiPts,'Curvature',[1,1],'EdgeColor','b');
    end
end
colormap(map);

set(ovFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(ovFig,strcat(ovPathname,ovFilename(1:end-4),'_Position'),'-dpdf'); 

stackOV = stack;

save(strcat(ovPathname,ovFilename(1:end-4),'.mat'),'stackOV'); 

%% Run the data for demonstrating excitation targeting
clc;
close all;
clear stack; clear stackXYfilt; clear stackXYfiltBGsub; clear ROIInfo; clear fullpath; clear positionDat;
[fullpath, stack, positionDat] = VRImDatLoad;
stackXYfilt = GaussFilt(stack);
stackXYfiltBGsub = BGSub(stackXYfilt);
ROIInfo = ROIExtract(stackXYfiltBGsub);

% Calculate the rotational and forward velocities
tStep = mean(diff(positionDat.t));
DeltaRot = diff(positionDat.OffsetRot);
notHuge = find(abs(DeltaRot) < 10);
OffsetRotDPlot = abs(smooth(DeltaRot(notHuge),100))./tStep;

OffsetForSmooth = smooth(positionDat.OffsetFor,100);
OffsetLatSmooth = smooth(positionDat.OffsetLat,100);
netSpeed = sqrt(diff(OffsetForSmooth).^2+diff(OffsetLatSmooth).^2)./tStep;

% Plot relevant data
cc=hsv(length(ROIInfo)+1);
cc(1,:) = 0;
allROIs = zeros(size(ROIInfo(1).logicalROI));

ROIFig = figure('units','normalized','outerposition',[0 0 1 1])
for roiplots = 1:length(ROIInfo)
    subplot(length(ROIInfo)+2,4,[1:3]+4*(roiplots));
    plot(positionDat.t(1)-(positionDat.tVR(1)-positionDat.tFrameGrab)/10000,ROIInfo(roiplots).ROIave,'color',cc(roiplots+1,:));
    hold on;
    scatter(positionDat.t(1)-(positionDat.tVR(1)-positionDat.tStim)/10000,zeros(1,length(positionDat.tStim)),'r')
    axis tight;
    axis off;
    xlim([positionDat.t(1) positionDat.t(end)]);
    allROIs = allROIs + ROIInfo(roiplots).logicalROI.*roiplots;
end
subplot(length(ROIInfo)+2,4,[4 8 12]);
imshow(allROIs);
axis equal;
axis off;
colormap(cc);
caxis([-0.5 roiplots+0.5]);

subplot(length(ROIInfo)+2,4,[16 20 24]);
imshow(stackOV(:,:,1),[0 upperLim])
hold on;
for roiNow=5:5
    roiPts = ROIs{roiNow}.ROIs(1).drawings.position;
    if strcmp(ROIs{roiNow}.ROIs(1).drawings.class,'impoly');
        for roiPt = 1:length(roiPts)
            if roiPt == length(roiPts)
                line([roiPts(roiPt,1) roiPts(1,1)], [roiPts(roiPt,2) roiPts(1,2)],'Color',[1 0 0],'LineWidth',5);
            else
                line([roiPts(roiPt,1) roiPts(roiPt+1,1)], [roiPts(roiPt,2) roiPts(roiPt+1,2)],'Color',[1 0 0],'LineWidth',5);
            end
        end
    else
        rectangle('Position',roiPts,'Curvature',[1,1],'EdgeColor','r','LineWidth',5);
    end
end
colormap(cc);

subplot(length(ROIInfo)+2,4,[1:3]);
plot(positionDat.t(1:end-1),netSpeed,'k');
axis tight;
set(gca,'FontSize',16);
ylabel('V_{forward} (cm/sec)','FontSize',16);

subplot(length(ROIInfo)+2,4,[5:7]+4*(length(ROIInfo)));
plot(positionDat.t(notHuge),OffsetRotDPlot,'k');
axis tight;
set(gca,'FontSize',16);
xlabel('Time (sec)','FontSize',16);
ylabel('V_{rot} (degrees/sec)','FontSize',16);

%% Look at the forward velocity alone after triggering
VfFig = figure('units','normalized','outerposition',[0 0 1 1])
plot(positionDat.t(1:end-1),netSpeed,'k');
hold on;
scatter(positionDat.t(1)-(positionDat.tVR(1)-positionDat.tStim)/10000,zeros(1,length(positionDat.tStim)),'r')
axis tight;
set(gca,'FontSize',16);
ylabel('V_{forward} (cm/sec)','FontSize',16);
xlabel('Time (sec)');

%% Save the data
set(ROIFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(ROIFig,strcat(fullpath{1}(1:end-4)),'-dpdf'); 

set(VfFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(VfFig,strcat(fullpath{1}(1:end-4),'_Vf'),'-dpdf');


save(strcat(fullpath{1}(1:end-4),'.mat'),'stackXYfiltBGsub','positionDat','fullpath','ROIInfo');

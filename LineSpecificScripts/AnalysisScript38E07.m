%% Run the data for demonstrating excitation targeting
clc;
close all;
clear;
[fullpath, stackMean, stackMaxInt, positionDat] = VRImDatLoad;
% stackXYfilt = GaussFilt(stack);
stackBGsub = BGSub(stackMaxInt);
% ROIInfo = ROIExtract(stackXYfiltBGsub);
ROIInfo = ROIExtract(stackBGsub);

% Calculate the rotational and forward velocities
[posRot, posFor, posLat] = PositionConverter(positionDat);
tStep = mean(diff(positionDat.t));
% DeltaRot = diff(positionDat.OffsetRot);
DeltaRot = diff(posRot);
notHuge = find(abs(DeltaRot) < 10);
OffsetRotDPlot = abs(smooth(DeltaRot(notHuge),100))./tStep;

% OffsetForSmooth = smooth(positionDat.OffsetFor,100);
% OffsetLatSmooth = smooth(positionDat.OffsetLat,100);
OffsetForSmooth = smooth(posFor,100);
OffsetLatSmooth = smooth(posLat,100);
netSpeed = sqrt(diff(OffsetForSmooth).^2+diff(OffsetLatSmooth).^2)./tStep;

% Plot relevant data
cc=hsv(length(ROIInfo)+1);
cc(1,:) = 0;
allROIs = zeros(size(ROIInfo(1).logicalROI));
num_planes = length(positionDat.tFrameGrab)/length(ROIInfo(1).ROIave);
ROIFig = figure('units','normalized','outerposition',[0 0 1 1])
for roiplots = 1:length(ROIInfo)
    subplot(length(ROIInfo)+2,4,[1:3]+4*(roiplots));
    plot(positionDat.t(1)-(positionDat.tVR(1)-positionDat.tFrameGrab(1:num_planes:end))/10000,ROIInfo(roiplots).ROIave,'color',cc(roiplots+1,:));
    hold on;
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

% subplot(length(ROIInfo)+2,4,[16 20 24]);
% imshow(mean(stack,3),[0 max(max(mean(stack,3)))]);
% colormap(cc);

subplot(length(ROIInfo)+2,4,[1:3]);
plot(positionDat.t(1:end-1),netSpeed,'k');
hold on;
plot(positionDat.t,positionDat.OffsetFor,'r');
axis tight;
set(gca,'FontSize',16);
ylabel('V_{forward} (cm/sec)','FontSize',16);

subplot(length(ROIInfo)+2,4,[5:7]+4*(length(ROIInfo)));
plot(positionDat.t(notHuge),OffsetRotDPlot,'k');
hold on;
plot(positionDat.t,positionDat.OffsetRot,'r');
axis tight;
set(gca,'FontSize',16);
xlabel('Time (sec)','FontSize',16);
ylabel('V_{rot} (degrees/sec)','FontSize',16);

 % Do correlation analysis
% figure;
% t1 = positionDat.tFrameGrab(1:num_planes:end)/10000;
% tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
% t2 = tSpan(1:length(OffsetRotDPlot));
% ROIMaxAct = max(ROIInfo.ROIave,[],1)./max(max(ROIInfo.ROIave));
% ROIforCC = interp1(t1,(ROIMaxAct-min(ROIMaxAct))/(max(ROIMaxAct)-min(ROIMaxAct)), t2);
% for i = 1:101
%     Acoef = 0.01*(i-1);
%     Bcoef = 0.01*(101-i);
%     testSig = Acoef*OffsetRotDPlot./max(OffsetRotDPlot)+Bcoef*smooth(netSpeed(1:length(OffsetRotDPlot)),360)./max(smooth(netSpeed,360));
%     checkVals = corrcoef(ROIforCC(~isnan(ROIforCC)),testSig(~isnan(ROIforCC)));
%     findMaxCC(i) = checkVals(2,1);
% end
% bestVal = find(findMaxCC == max(findMaxCC));
% Acoef = 0.01*(bestVal-1);
% Bcoef = 0.01*(101-bestVal);
% plot(positionDat.tFrameGrab(1:num_planes:end)/10000,(ROIMaxAct-min(ROIMaxAct))/(max(ROIMaxAct)-min(ROIMaxAct)),'r');
% axis tight;
% set(gca,'FontSize',16);
% hold on;
% plot(tSpan(1:length(OffsetRotDPlot)),Acoef*OffsetRotDPlot./max(OffsetRotDPlot)+Bcoef*smooth(netSpeed(1:length(OffsetRotDPlot)),360)./max(smooth(netSpeed,360)),'b');
% legend('Normalized Fmax',strcat(num2str(Acoef),'*V_{rot}+',num2str(Bcoef),'*V_{for}'));
% xlabel('Time');

%% Load the data
clc;
close all;
clear;
[fullpath, stackMean, stackMaxInt, positionDat] = VRImDatLoadDiscard;

%% Extract the ROIs
ROIInfo = ROIExtract(stackMaxInt);

%% Plot the ROIs over time as compared to the velocity and the 
% Normalize the ROIs
allROI = zeros(length(ROIInfo),length(ROIInfo(1).ROIave));
for ROIstep = 1:length(ROIInfo)
    ROIsort = sort(ROIInfo(ROIstep).ROIave);
    allROI(ROIstep,:) = (ROIInfo(ROIstep).ROIave-mean(ROIsort(1:end/10)))./(max(ROIInfo(ROIstep).ROIave)-mean(ROIsort(1:end/10)));
end

% Downsample the postions and find the velocity
tDS = zeros(1,length(ROIInfo(1).ROIave));
OffsetLatDS = zeros(1,length(ROIInfo(1).ROIave));
OffsetForDS = zeros(1,length(ROIInfo(1).ROIave));
OffsetRotDS = zeros(1,length(ROIInfo(1).ROIave));
for DSstep = 1:length(OffsetLatDS);
    tNow = positionDat.tFrameGrab(3+6*(DSstep-1))/10000-positionDat.tVR(1)/10000+positionDat.t(1);
    if tNow > positionDat.t(1)
        tDS(DSstep) = tNow;
        tSpan = find(positionDat.t>=tDS(DSstep));
        OffsetLatDS(DSstep) = positionDat.OffsetLat(tSpan(1));
        OffsetForDS(DSstep) = positionDat.OffsetFor(tSpan(1));
        OffsetRotDS(DSstep) = positionDat.OffsetRot(tSpan(1));
    end
end

vFor = smooth(abs(diff(sqrt(OffsetLatDS.^2+OffsetForDS.^2))./diff(tDS)),15);

% Find the PVA of the FB activity and the object location
angsraw = (1:(length(ROIInfo)-1))*2*pi/(length(ROIInfo)-1)-pi;
angsraw = angsraw';
for ts = 1:length(allROI)
    meanAngRaw(ts) = circ_mean(angsraw, squeeze(allROI(2:end,ts)));
    meanIntRaw(ts) = circ_r(angsraw, squeeze(allROI(2:end,ts)));
end

relAng = atan2(OffsetForDS,-OffsetLatDS);
netAngle(find(-OffsetLatDS > 0)) = -relAng(find(-OffsetLatDS > 0))-pi/2+pi/180*OffsetRotDS(find(-OffsetLatDS > 0));
netAngle(find(-OffsetLatDS < 0)) = -relAng(find(-OffsetLatDS < 0))+3*pi/2+pi/180*OffsetRotDS(find(-OffsetLatDS < 0));
netAngle=mod(netAngle+pi,2*pi)-pi;

% Find the minimum offset between the PVA and the object angle;
netAngleUnwrap = UnWrap(netAngle,2,0);

for offset = pi/32:pi/32:2*pi;
    netAngleOffset = mod(netAngle-offset,2*pi)-pi;
    angError(round(offset*32/pi)) = sum(min(abs(netAngleOffset-meanAngRaw),abs(2*pi-(netAngleOffset-meanAngRaw))));
end

minOffset = find(angError == min(angError));
minOffset = minOffset*pi/32;
netAngleOffset = mod(netAngle-minOffset,2*pi)-pi;

%% Plot interesting things

PVAvsFB =  figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,1,1)
plot(tDS(2:end),vFor,'g','LineWidth',2);
xlim([tDS(1) tDS(end)]);
ylim([0 1]);
set(gca,'FontSize',16);
ylabel('Speed (cm/s)')

subplot(4,1,2)
hold on;
plot(positionDat.tFrameGrab(3:6:end)/10000-positionDat.tVR(1)/10000+positionDat.t(1),allROI(1,:))
plot(positionDat.tFrameGrab(3:6:end)/10000-positionDat.tVR(1)/10000+positionDat.t(1),mean(allROI(2:end,:),1),'k');
xlim([tDS(1) tDS(end)]);
ylim([0 1]);
set(gca,'FontSize',16);
ylabel('FB Intensity (\Delta F/F)');

subplot(4,1,3)
hold on;
contourf(positionDat.tFrameGrab(3:6:end)/10000-positionDat.tVR(1)/10000+positionDat.t(1),[-pi angsraw'],vertcat(allROI(end,:),allROI(2:end,:)),'LineStyle','None')
% plot(tDS,smooth(meanAngRaw,10),'k');
set(gca,'FontSize',16);

subplot(4,1,4);
hold on;
plot(tDS,smooth(meanAngRaw,10),'k');
plot(tDS,netAngleOffset,'g');
% plot(tDS,OffsetRotDS./180*pi,'r');
xlim([tDS(1) tDS(end)]);
ylim([-pi pi]);
set(gca,'FontSize',16);
ylabel('Predicted PVA (rad)');
xlabel('Time (sec)');

set(PVAvsFB,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PVAvsFB,strcat(fullpath{1}(1:end-4),'_Comparisons'),'-dpdf'); 

%% Tie in the position
figure('units','normalized','outerposition',[0 0 1 1]);
colorsDist = vertcat([1 0.5 0], [1 0 0], [1 1 0], [0 1 0], [0 1 1], [0 0 1], [1 0 1], [1 1 1]);
for ROIit = 2:9
    cNow = vertcat((1-squeeze(allROI(ROIit,:))).*colorsDist(ROIit-1,1), (1-squeeze(allROI(ROIit,:))).*colorsDist(ROIit-1,2), (1-squeeze(allROI(ROIit,:))).*colorsDist(ROIit-1,3));
    subplot(3,3,ROIit)
    scatter(OffsetLatDS,OffsetForDS,4,cNow');
    viscircles([0 0], 1,'EdgeColor','b');
    axis equal;
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,[1:2]);
hold on;
for pt = 1:length(OffsetLatDS)
    ROIit = 1+find(allROI(2:end,pt) == max(allROI(2:end,pt)));
    scatter(OffsetLatDS(pt),OffsetForDS(pt),10,colorsDist(ROIit-1,:),'filled')
end
viscircles([0 0], 1,'EdgeColor','b');
axis equal;

subplot(1,3,3);
scatter([0 0 0 0 0 0 0 0],[1:8],50,colorsDist,'filled');
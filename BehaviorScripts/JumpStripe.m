%% Change to the appropriate directory
allPathname = uigetdir('D:/Imaging','Select the directory');
cd(allPathname);

%% Plot things!
clear;
close all;
clc;
flyFilename = uigetfile;
load(flyFilename);

num_ROIs = size(ROIaveMax,1);
num_planes = round(length(positionDat.tFrameGrab)/length(ROIaveMax));

% Pull out the trial periods
trialbreaks = find(abs(diff(positionDat.direction(1:end-2)))>=1);
trialbreaks = vertcat(0,trialbreaks,length(positionDat.t));
numTrials = length(trialbreaks)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the circular mean
angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
angsraw = angsraw';
clear meanAngRaw;
clear meanIntRaw;
for ts = 1:length(ROIaveMax)
    meanAngRaw(ts) = circ_mean(angsraw, squeeze(ROIaveMax(:,ts)));
    meanIntRaw(ts) = circ_r(angsraw, squeeze(ROIaveMax(:,ts)));
end

tSpan = positionDat.t-positionDat.t(1)+positionDat.tVR(1)/10000;
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);

% Downsample the data to 20Hz for consistency with Hannah
sampFreq = 1/mean(diff(positionDat.t));
pts2sum = round(sampFreq/20);
tDS = zeros(floor(length(positionDat.OffsetLat)/pts2sum),1);
OffsetRotDS = tDS;
OffsetForDS = tDS;
OffsetLatDS = tDS;
for DSstep = 1:length(tDS);
    tDS(DSstep) = mean(positionDat.t(1+pts2sum*(DSstep-1):pts2sum*DSstep));
    OffsetRotDS(DSstep) = mean(positionDat.OffsetRot(1+pts2sum*(DSstep-1):pts2sum*DSstep));
    OffsetForDS(DSstep) = mean(positionDat.OffsetFor(1+pts2sum*(DSstep-1):pts2sum*DSstep));
    OffsetLatDS(DSstep) = mean(positionDat.OffsetLat(1+pts2sum*(DSstep-1):pts2sum*DSstep));
end

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFG-minFG+1,2);
netAngleMatch = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
end

OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
% Set the timestamps for the PVA
PVASpan = OffsetRotMatch(:,1)-positionDat.t(1)+positionDat.tVR(1)/10000;

flyOverview = figure('units','normalized','outerposition',[0 0 1 1]);
cylDiameter = 1;
subplot(5,3,[1 4]);
hold on;
scatter(-OffsetLatDS,OffsetForDS,1);
for plotPt = 1:length(tDS)
   line([-OffsetLatDS(plotPt) -OffsetLatDS(plotPt) - 0.15*sin(OffsetRotDS(plotPt)*pi/180)],[OffsetForDS(plotPt) OffsetForDS(plotPt) + 0.15*cos(OffsetRotDS(plotPt)*pi/180)]); 
end
axis equal;
scatter(0,-7,100,'p','k')
xlabel('Distance (cm)');
ylabel('Distance (cm)');
set(gca,'FontSize',16);

% ***NOTE*** The timescale of the DAQ data hs been chosen rather than the
% program output timestamps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the activity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(5,3,2:3);
hold on;
imagesc('XData',positionDat.tFrameGrab(1:num_planes:end)/10000,'YData',2*pi*((1:num_ROIs)-(num_ROIs+1)/2)/(num_ROIs-1),'CData',ROIaveMax(:,1:ceil(length(positionDat.tFrameGrab)/num_planes)));
colormap(brewermap(64, 'Blues'));
plot(positionDat.tFrameGrab(1:num_planes:end)/10000,meanAngRaw(:,1:ceil(length(positionDat.tFrameGrab)/num_planes)),'k');
caxis([min(min(ROIaveMax)) max(max(ROIaveMax))]);
set(gca,'FontSize',16);
ylabel('Position (rad)');
title('Maximum Intensity');
ylim([-pi pi]);
xlim([PVASpan(1) PVASpan(end)])
for trial = 2:numTrials
    line([positionDat.t(trialbreaks(trial))-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000 positionDat.t(trialbreaks(trial))-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000], [-pi pi], 'LineStyle','--');
end

% Plot the PVA amplitude
% Code the PVA amplitude by color
PVAAmpCol = (1-meanIntRaw(:,1:ceil(length(positionDat.tFrameGrab)/num_planes))./max(meanIntRaw(:,1:ceil(length(positionDat.tFrameGrab)/num_planes))))'*[1 1 1];
scatter(positionDat.tFrameGrab(1:num_planes:end)/10000,zeros(length(meanIntRaw(:,1:ceil(length(positionDat.tFrameGrab)/num_planes))),1)-pi,20,PVAAmpCol,'filled');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the circular mean with the global rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanAngRawUnwrap = UnWrap(meanAngRaw(minFG:maxFG),1.5,1);

subplot(5,3,5:6);
hold on;
plot(OffsetRotMatch(:,1)-positionDat.t(1)+positionDat.tVR(1)/10000,OffsetRotMatchUnwrap - mean(OffsetRotMatchUnwrap),'--g');
plot(OffsetRotMatch(:,1)-positionDat.t(1)+positionDat.tVR(1)/10000,OffsetRotMatchUnwrap*1.5 - mean(OffsetRotMatchUnwrap*1.5),'g');
plot(PVASpan,meanAngRawUnwrap - mean(meanAngRawUnwrap),'k');
xlim([PVASpan(1) PVASpan(end)]);
ylabel('Body Rotation (rad)','FontSize',16);
set(gca,'FontSize',16);

for trial = 2:numTrials
    line([positionDat.t(trialbreaks(trial))- positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(trialbreaks(trial))- positionDat.t(1)+positionDat.tVR(1)/10000], [-pi pi], 'LineStyle','--');
end

startCorr = find(OffsetRotMatch(:,1) > 15 + positionDat.t(trialbreaks(2)) - positionDat.t(1)+positionDat.tVR(1)/10000);
stopCorr = find(OffsetRotMatch(:,1) < positionDat.t(trialbreaks(3)) - positionDat.t(1)+positionDat.tVR(1)/10000);
CCs = corrcoef(OffsetRotMatchUnwrap(startCorr(1):stopCorr(end)),meanAngRawUnwrap(startCorr(1):stopCorr(end)));
text(positionDat.t(trialbreaks(2))- positionDat.t(1)+positionDat.tVR(1)/10000, pi, strcat('Correlation = ', num2str(CCs(1,2))),'FontSize',16); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the running difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,8:9);
% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 5; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(abs(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt))));
end

hold on;
plot(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,'k');
scatter(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,20,PVAAmpCol(minFG:length(runDiffGlob)+minFG-1,:),'filled');
set(gca,'FontSize',14);
title(['Running Difference (',num2str(tDiff), ' sec)']);
ylabel('Radians');
axis tight;
xlim([PVASpan(1) PVASpan(end)])
ylim([0 pi])

subplot(5,3,11:12);
% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 10; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(abs(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt))));
end

hold on;
plot(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,'k');
scatter(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,20,PVAAmpCol(minFG:length(runDiffGlob)+minFG-1,:),'filled');
set(gca,'FontSize',14);
title(['Running Difference (',num2str(tDiff), ' sec)']);
ylabel('Radians');
axis tight;
xlim([PVASpan(1) PVASpan(end)])
ylim([0 pi])

subplot(5,3,14:15);
% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 20; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(abs(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt))));
end

hold on;
plot(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,'k');
scatter(OffsetRotMatch(1:end-nPtsDiff+1,1),runDiffGlob,20,PVAAmpCol(minFG:length(runDiffGlob)+minFG-1,:),'filled');
set(gca,'FontSize',14);
title(['Running Difference (',num2str(tDiff), ' sec)']);
ylabel('Radians');
axis tight;
xlim([PVASpan(1) PVASpan(end)])
ylim([0 pi])

%% Plot the periods where the stripe is invisible
angShift = 7.5 * pi/180;
darkPeriods = find(abs(OffsetRotMatch(:,2)+angShift) > 2*pi/3+angShift);
subplot(5,3,11:12);
scatter(OffsetRotMatch(darkPeriods,1)-OffsetRotMatch(1,1)+positionDat.tVR(1)/10000,3+zeros(length(darkPeriods),1),'k');

%% Plot the periods where the stripe jumps
bumpJumps = find(abs(diff(OffsetRotMatch(:,2)))>pi/2);
subplot(5,3,11:12);
scatter(OffsetRotMatch(bumpJumps,1),3+zeros(length(bumpJumps),1),'k','filled');

%% Save the plot
print(flyOverview,strcat(flyFilename(1:end-4),'_Overview'),'-djpeg'); 
set(flyOverview,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(flyOverview,strcat(flyFilename(1:end-4),'_Overview'),'-dpdf'); 

%% Plot some stats
runDiffComp = figure('units','normalized','outerposition',[0 0 1 1]);

% Find the net difference between the stripe and the normal and scaled
% versions
% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 10; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(abs(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt))));
end

runDiffGlobSigned = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlobSigned(runIt) = mean(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt)));
end

runDiffGlobScaled = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlobScaled(runIt) = mean(1.5*OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-1.5*OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt)));
end

subplot(2,2,1);
text(0,1,'Stats:','FontSize',16);
text(0,0.5,strcat('Sum of diff: ',num2str(sum(runDiffGlobSigned))),'FontSize',16);
text(0,0,strcat('Sum of scaled diff: ',num2str(sum(runDiffGlobScaled))),'FontSize',16);
axis off;



startFit = find(OffsetRotMatch(:,1)-OffsetRotMatch(1,1) > 15 + positionDat.t(trialbreaks(2))-positionDat.t(1));
stopFit = find(OffsetRotMatch(:,1)-OffsetRotMatch(1,1) < positionDat.t(trialbreaks(3))-positionDat.t(1));
startFitPVA = find(positionDat.tFrameGrab/10000 > 15 + positionDat.t(trialbreaks(2))-positionDat.t(1)+positionDat.tVR(1)/10000);
numFramesFit = stopFit(end)-startFit(1);

% Only consider 15 sec in to the trial until the trial end
calcSpan = ceil(startFitPVA(1)/num_planes):min(ceil(startFitPVA(1)/num_planes)+numFramesFit,length(runDiffGlob));

OffsetRotMatch = zeros(maxFG-minFG+1,2);
for interp = minFG:maxFG
    tMatch = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
end

% Find mean values
nBinsPos = 37;
posAngles = linspace(-pi,pi,nBinsPos);
posValsRot = zeros(size(posAngles));
posBinsRot = zeros(size(posAngles));
for posIt = 1:min(numFramesFit,length(runDiffGlob)-ceil(startFitPVA(1)/num_planes)+1)
    if meanIntRaw(ceil(startFitPVA(1)/num_planes)+posIt-1)./max(meanIntRaw) > 0.5 % Ignore PVAs that are less than half the max
        spotRot = ceil(OffsetRotMatch(ceil(startFitPVA(1)/num_planes)+posIt-1,2)*(nBinsPos)/(2*pi)+(nBinsPos)/2);
        posValsRot(spotRot) = posValsRot(spotRot) + runDiffGlob(ceil(startFitPVA(1)/num_planes)+posIt-1);
        posBinsRot(spotRot) = posBinsRot(spotRot) + 1;
    end
end

% Scatter plots and means;
subplot(2,2,2);
hold on;
scatter(OffsetRotMatch(calcSpan,2), runDiffGlob(calcSpan), 10, PVAAmpCol(calcSpan,:),'filled');
xlabel('Global rotation (rad)');
ylabel('Running difference');
set(gca,'FontSize',16);
plot(posAngles,posValsRot./posBinsRot,'k','LineWidth',2);
xlim([-pi pi]);
title('Global rotation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at the walking S.D. vs. correlation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runCorr = zeros(length(meanAngRawUnwrap)-199,1);
runSD = zeros(length(meanAngRawUnwrap)-199,1);
for slide = 1:25:length(runCorr)
    slideCorr = corrcoef(OffsetRotMatchUnwrap(slide:slide+199),meanAngRawUnwrap(slide:slide+199));
    runCorr(slide) = slideCorr(1,2);
    runSD(slide) = std(OffsetRotMatchUnwrap(slide:slide+199));
end

corrSD = hist2d([runSD, runCorr], linspace(0, 3.5, 35), linspace(-1, 1, 20));

subplot(2,2,4);
pc = pcolor(linspace(-1, 1, 19), linspace(0, 3.5, 34),corrSD);
set(pc, 'EdgeColor', 'none');
colormap(map);
axis equal;
axis tight;
xlabel('Correlation');
ylabel('Walking rotation S.D. (rad)');
set(gca,'FontSize',16); 

%% Save the plot
print(runDiffComp,strcat(flyFilename(1:end-4),'_DiffStats'),'-djpeg'); 
set(flyOverview,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(runDiffComp,strcat(flyFilename(1:end-4),'_DiffStats'),'-dpdf'); 

%% Look at error around jumps
%%%%%%%%% Need to account for min FG
tBuffer = 10;
tAdd = round(tBuffer/(2*tStep));
jumpTimes = [];
jumpIt = 1;

for i=startCorr(1):stopCorr(end)
    for jump = 1:length(bumpJumps)
       if bumpJumps(jump)-tAdd <= i & i <= bumpJumps(jump)+tAdd
           jumpTimes(jumpIt) = i;
           jumpIt = jumpIt + 1;
       end
    end
end
jumpTimes = unique(jumpTimes);
otherTimes = setdiff([startCorr(1):stopCorr(end)],jumpTimes);

% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 5; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt)));
end

histBins = linspace(-1.5,1.5,41);

errorCmp = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
[nb,xb] = hist(runDiffGlob(jumpTimes),histBins);
bh = bar(histBins,nb);
set(gca,'FontSize',16);
xlim([-2 2]);
title('Jump periods');

subplot(2,1,2);
[nbOther,xbOther] = hist(runDiffGlob(otherTimes),histBins);
bh = bar(histBins,nbOther);
set(gca,'FontSize',16);
xlim([-2 2]);
xlabel('Running difference (rad)');
title('Other times');

nbOff = [];
nbOn = [];

%% Look at error in and around dark periods
tBuffer = 10;
darkPeriods(end) = length(OffsetRotMatch);
tAdd = round(tBuffer/(2*tStep));
darkTimes = [];
darkOnTimes = [];
darkOffTime = [];
darkIt = 1;
darkOn = find(diff(diff(darkPeriods()))<0)+1;
darkOnIt = 1;
darkOff = find(diff(diff(darkPeriods))>0)+1;
darkOffIt = 1;

for i=startCorr(1)+tAdd:stopCorr(end)
    for dark = 1:length(darkPeriods)
       if darkPeriods(dark)-tAdd <= i & i <= darkPeriods(dark)+tAdd
           darkTimes(darkIt) = i;
           darkIt = darkIt + 1;
       end
    end
    for darkOnStep = 1:length(darkOn)
       if darkPeriods(darkOn(darkOnStep))-tAdd <= i & i <= darkPeriods(darkOn(darkOnStep))+tAdd
           darkOnTimes(darkOnIt) = i;
           darkOnIt = darkOnIt + 1;
       end
    end
    for darkOffStep = 1:length(darkOff)
       if darkPeriods(darkOff(darkOffStep))-tAdd <= i & i <= darkPeriods(darkOff(darkOffStep))+tAdd
           darkOffTimes(darkOffIt) = i;
           darkOffIt = darkOffIt + 1;
       end
    end
end
darkTimes = unique(darkTimes);
otherTimes = setdiff([startCorr(1):stopCorr(end)],darkTimes);
darkOnTimes = unique(darkOnTimes);
otherOnTimes = setdiff([startCorr(1):stopCorr(end)],darkOnTimes);
darkOffTimes = unique(darkOffTimes);
otherOffTimes = setdiff([startCorr(1):stopCorr(end)],darkOffTimes);

figure;
hold on;
scatter(OffsetRotMatch(darkPeriods,1),zeros(length(darkPeriods),1),'k');
scatter(OffsetRotMatch(darkTimes,1),zeros(length(darkTimes),1)+1,'k');
scatter(OffsetRotMatch(darkPeriods(darkOn),1),zeros(length(darkOn),1),'g','filled');
scatter(OffsetRotMatch(darkOnTimes,1),zeros(length(darkOnTimes),1)+2,'g','filled');
scatter(OffsetRotMatch(darkPeriods(darkOff),1),zeros(length(darkOff),1),'r','filled');
scatter(OffsetRotMatch(darkOffTimes,1),zeros(length(darkOffTimes),1)+3,'r','filled');


% Find the running difference
tStep = mean(diff(OffsetRotMatch(:,1)));
tDiff = 5; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotMatch)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotMatch)-nPtsDiff+1
    runDiffGlob(runIt) = mean(OffsetRotMatchUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotMatchUnwrap(runIt)-(meanAngRawUnwrap(runIt:runIt+nPtsDiff-1)-meanAngRawUnwrap(runIt)));
end

% Set up the bins
histBins = linspace(-1.5,1.5,41);

errorCmp = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
[nb,xb] = hist(runDiffGlob(darkTimes),histBins);
bh = bar(histBins,nb);
set(bh,'facecolor',[0.3 0.3 0.3]);
hold on;

[nbOff,xbOff] = hist(runDiffGlob(darkOffTimes),histBins);
bhOff = bar(histBins,nbOff);
set(bhOff,'facecolor',[1 0 0]);
alpha(bhOff,0.5);

[nbOn,xbOn] = hist(runDiffGlob(darkOnTimes),histBins);
bhOn = bar(histBins,nbOn);
set(bhOn,'facecolor',[0 1 0]);
alpha(bhOn,0.5);

set(gca,'FontSize',16);
xlim([-1.5 1.5]);
title('Dark periods');


subplot(2,1,2);
[nbOther,xbOther] = hist(runDiffGlob(otherTimes),histBins);
bhOther = bar(histBins,nbOther);
set(bhOther,'facecolor',[0.5 0.5 1]);
set(gca,'FontSize',16);
xlim([-1.5 1.5]);
xlabel('Running difference (rad)');
title('Other times');

%% Save the plot and the data
print(errorCmp,strcat(flyFilename(1:end-4),'_ErrorCmp'),'-djpeg'); 
set(flyOverview,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(errorCmp,strcat(flyFilename(1:end-4),'_ErrorCmp'),'-dpdf'); 

save(strcat(flyFilename(1:end-4),'_ErrorCmp'),'histBins','nb','nbOff','nbOn','nbOther');
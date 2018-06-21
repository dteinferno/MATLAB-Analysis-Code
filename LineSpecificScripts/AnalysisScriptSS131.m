%% Clear out old data and load the fresh stuff
clc;
close all;
clear;
[fullpath, stackMean, stackMaxInt, positionDat] = VRImDatLoad;
% [fullpath, stackMean, stackMaxInt, positionDat] = VRImDatLoadDiscard;
% stackXYfilt = GaussFilt(stack);
% stackXYfiltBGsub = BGSub(stackXYfilt);
% sgolayOrder = 3;
% sgolayWindow = 7;
% stackXYTfiltBGsub = sgolayfilt(stackXYfiltBGsub,sgolayOrder,sgolayWindow,[],3);

%% Extract the ROIs and find the DF/F
% Fit an ellipse over the entire EB and then break it into ROI sections
A = mean(stackMaxInt,3);
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

ROIs = zeros(num_ROIs,256,256);
for incROI = 1:num_ROIs
    xROI(1) = ellipseCenter(1);
    yROI(1) = ellipseCenter(2);
    for angs=1:5
        xROI(angs+1) = ellipseCenter(1)+ellipseWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        yROI(angs+1) = ellipseCenter(2)+ellipseHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
    end
    ROIs(incROI,:,:) = roipoly(A, xROI,yROI);
end

% Plot the ROI on each figure and calculate the average for both the
% Maximum Intesity and Mean projections
ROIaveREFMax = zeros(num_ROIs,length(stackMaxInt));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:length(stackMaxInt)
    if mod(i,100)==0
        waitbar(i/length(stackMaxInt),h,['Calculating frame# ' num2str(i) ' out of ' num2str(length(stackMaxInt))]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(stackMaxInt(:,:,i));
        ROI = A(logical(squeeze(ROIs(incROI,:,:))));
        ROIaveREFMax(incROI,i) = mean2(ROI);
    end
end
delete(h);
ROIaveMax = zeros(num_ROIs,length(stackMaxInt));
for incROI = 1:num_ROIs
    ROILow = sort(squeeze(ROIaveREFMax(incROI,:)));
    ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./sum(ROILow(1:floor(end/10)));
end

ROIaveREFMean = zeros(num_ROIs,length(stackMean));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:length(stackMean)
    if mod(i,100)==0
        waitbar(i/length(stackMean),h,['Calculating frame# ' num2str(i) ' out of ' num2str(length(stackMean))]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(stackMean(:,:,i));
        ROI = A(logical(squeeze(ROIs(incROI,:,:))));
        ROIaveREFMean(incROI,i) = mean2(ROI);
    end
end
delete(h);
ROIaveMean = zeros(num_ROIs,length(stackMean));
for incROI = 1:num_ROIs
    ROILow = sort(squeeze(ROIaveREFMean(incROI,:)));
    ROIaveMean(incROI,:) = ROIaveREFMean(incROI,:)./sum(ROILow(1:floor(end/10)));
end

clear ellipse* xROI yROI

%% Calculate the circular mean and plot it on top of the activity
num_planes = ceil(length(positionDat.tFrameGrab)/length(ROIaveMax));
% Discriminate between ther different trials
trialbreaks1 = find(abs(diff(positionDat.closed))>=1);
trialbreaks2 = find(abs(diff(positionDat.direction))>=1);
trialbreaks = union(trialbreaks1,trialbreaks2);
trialbreaks = vertcat(1,trialbreaks,length(positionDat.t));
numTrials = length(trialbreaks)-1;
clear trialbreaks1 trialbreaks2

act = figure('units','normalized','outerposition',[0 0 1 1]);

for projs = 1:2
    if projs == 1
        ROIave = ROIaveMax;
    else
        ROIave = ROIaveMean;
    end
    % Calculate and plot the circular mean
    angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
    angsraw = angsraw';
    clear meanAngRaw;
    clear meanIntRaw;
    for ts = 1:length(ROIave)
        meanAngRaw(ts) = circ_mean(angsraw, squeeze(ROIave(:,ts)));
        meanIntRaw(ts) = circ_r(angsraw, squeeze(ROIave(:,ts)));
    end

    meanAngRawPlot = meanAngRaw;
    for i=1:length(meanAngRaw)-1
        if abs(meanAngRaw(i)-meanAngRaw(i+1)) > pi
            meanAngRawPlot(i) = NaN;
        end
    end

    % Plot the activity and circular mean
    subplot(3,4,(1:4)+4*(projs-1));
    contourf(positionDat.tFrameGrab(1:num_planes:end)/10000,2*pi*((1:num_ROIs)-(num_ROIs+1)/2)/(num_ROIs-1),ROIave(:,1:ceil(length(positionDat.tFrameGrab)/num_planes)),'LineStyle','none');
    map = cbrewer('seq','Blues',64);
    colormap(map);
    caxis([min(min(ROIave)) max(max(ROIave))]);
    hold on;
    plot(positionDat.tFrameGrab(1:num_planes:end)/10000,meanAngRawPlot(1:ceil(length(positionDat.tFrameGrab)/num_planes)),'k');
    set(gca,'FontSize',16);
    ylabel('Position (rad)');
    if projs == 1
        title('Maximum Intesity');
    else
        title('Mean');
        xlabel('Time (sec)');
    end 
    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
end

maxORmean = input('Use max [1] or mean [0]?');
if maxORmean
    ROIave = ROIaveMax;  
    clear meanAngRaw;
    clear meanIntRaw;
    for ts = 1:length(ROIave)
        meanAngRaw(ts) = circ_mean(angsraw, squeeze(ROIave(:,ts)));
        meanIntRaw(ts) = circ_r(angsraw, squeeze(ROIave(:,ts)));
    end
end

clear A angs 

%% Unwrap the signals for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prep the data for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smooth the data
OffsetRotPlot =pi/180* smooth(positionDat.OffsetRot,10);

% Find the object orientation
relAng = atan2(positionDat.OffsetFor,-positionDat.OffsetLat);
netAngle(find(-positionDat.OffsetLat > 0)) = -relAng(find(-positionDat.OffsetLat > 0))-pi/2+pi/180*positionDat.OffsetRot(find(-positionDat.OffsetLat > 0));
netAngle(find(-positionDat.OffsetLat < 0)) = -relAng(find(-positionDat.OffsetLat < 0))+3*pi/2+pi/180*positionDat.OffsetRot(find(-positionDat.OffsetLat < 0));
netAngle=mod(netAngle+pi,2*pi)-pi;


% Downsample the rotation to match the imaging data
minFG = ceil(max(find(positionDat.tFrameGrab <= positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);
OffsetRotDS = zeros(maxFG-minFG+1,2);
netAngleDS = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tDS = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    OffsetRotDS(interp-minFG+1,1) = positionDat.t(tDS(1));
    OffsetRotDS(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tDS(1));
    netAngleDS(interp-minFG+1) = netAngle(tDS(1));
end
netAngle = smooth(netAngle,10);
    
% Unwrap the rotation
OffsetRotDSUnwrap = UnWrap(OffsetRotDS(:,2),2);

% Unwrap the object orientation
netAngleDSUnwrap = UnWrap(netAngleDS,2);

% Unwrap the activity
meanAngBodUnwrap = UnWrap(meanAngRaw(minFG:maxFG),1.5);

%% Compare the PV and the movement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the PV to the angle to the fly's movements and environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotComp = figure('units','normalized','outerposition',[0 0 1 1]);
c = linspace(1, 10, length(positionDat.OffsetLat(1:10:end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the fly's position over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,4,1);
hold on;
scatter(positionDat.OffsetLat(1:10:end),positionDat.OffsetFor(1:10:end),1,c);
if (strcmp(positionDat.exType,'LightCylFlatBack') || strcmp(positionDat.exType,'LightCyl3ObjBack'))
    viscircles([0 0], 1,'EdgeColor','b');
end
if (strcmp(positionDat.exType,'NoCyl3ObjBack') || strcmp(positionDat.exType,'LightCyl3ObjBack'))
    line([0 2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 0.5 0.5]);
    line([0 -2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 1 0]);
    line([0 -2*cos(pi/3)],[0 -2*sin(pi/3)],'LineStyle','--','Color',[1 0 0]);
end
axis equal;
set(gca,'FontSize',16);
xlabel('cm');ylabel('cm');


if strcmp(positionDat.exType,'NoCylFlatBack') || strcmp(positionDat.exType,'Collada file used: D:\Environments\SineBG')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the rotation to the activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    subplot(3,4,5:12);
    hold on;
    tprog = zeros(1,length(c))+min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap));
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1),OffsetRotDSUnwrap,'g');
    plot(OffsetRotDS(:,1),meanAngBodUnwrap,'k','LineWidth',1.5);
    axis tight;
    ylabel('Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap)) max(max(OffsetRotDSUnwrap),max(meanAngBodUnwrap))], 'LineStyle','--');
    end
    colormap('parula');
    
    for trial = 1:numTrials-1
        tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
        cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
        text(OffsetRotDS(tLimits(1),1),mean(OffsetRotDSUnwrap), ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
        set(gca,'FontSize',16); 
    end   
end

if strcmp(positionDat.exType,'LightCylFlatBack')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shift data for plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PVAOffsetBodAll = meanAngRaw(minFG:maxFG) - OffsetRotDS(:,2)';
    PVAOffsetObjAll = meanAngRaw(minFG:maxFG) - netAngleDS';
    for it = 1:length(meanAngRaw(minFG:maxFG))
        if PVAOffsetBodAll(it) < 0
            PVAOffsetBodAll(it) = PVAOffsetBodAll(it) + 2*pi;
        end
        if PVAOffsetObjAll(it) < 0
            PVAOffsetObjAll(it) = PVAOffsetObjAll(it) + 2*pi;
        end
    end
    PVAOffsetBodAll(isnan(PVAOffsetBodAll)) = 0;
    PVAOffsetObjAll(isnan(PVAOffsetObjAll)) = 0;
    PVAOffsetBod = mean(PVAOffsetBodAll);
    PVAOffsetObj = mean(PVAOffsetObjAll);

    % Shift the PVA by the appropriate amount 
    meanAngBod = meanAngRaw(minFG:maxFG) - PVAOffsetBod;
    meanAngBod = mod(meanAngBod-pi,2*pi)-pi;
    meanAngObj = meanAngRaw(minFG:maxFG) - PVAOffsetObj;
    meanAngObj = mod(meanAngObj-pi,2*pi)-pi;
    
    PVASpan = OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    
    % Get rid of wrap-around artifacts
    for i=1:length(OffsetRotDS)-1
         if abs(OffsetRotDS(i,2)-OffsetRotDS(i+1,2)) > pi/2
            OffsetRotDS(i,2) = NaN;
         end
         if abs(netAngle(i)-netAngle(i+1)) > pi/2
            netAngle(i) = NaN;
         end
         if abs(meanAngBod(i)-meanAngBod(i+1)) > pi/2
            meanAngBod(i) = NaN;
         end
         if abs(meanAngObj(i)-meanAngObj(i+1)) > pi/2
            meanAngObj(i) = NaN;
         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the local and global rotation to the activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    subplot(3,4,5:8);
    hold on;
    tprog = zeros(1,length(c))-pi;
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000,OffsetRotDS(:,2),'g');
    plot(PVASpan,meanAngBod,'k');
    axis tight;
    ylabel('Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(tSpan(1),0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
    
    subplot(3,4,9:12);
    hold off;
    CylWidth = 2*atan(1./sqrt(positionDat.OffsetFor.^2+positionDat.OffsetLat.^2));
    errorbarPatch2(tSpan,netAngle, CylWidth,[0 0 1]);
    plot(PVASpan,meanAngObj,'k');
    axis tight;
    ylabel('Object Orientation (rad)','FontSize',16);
    xlabel('Time (sec)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    colormap('parula');
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(netAngleDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(positionDat.tVR(1)/10000,0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
    p1=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[2*pi/3 pi pi 2*pi/3],'k');
    p2=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[-2*pi/3 -pi -pi -2*pi/3],'k');
    set(p1,'FaceAlpha',0.25); 
    set(p2,'FaceAlpha',0.25); 
    ylim([-pi pi]);
    
    
    
end

if strcmp(positionDat.exType,'NoCylFlatBack')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the rotation to the activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    subplot(3,4,5:12);
    hold on;
    tprog = zeros(1,length(c))+min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap));
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1),OffsetRotDSUnwrap,'g');
    plot(OffsetRotDS(:,1),meanAngBodUnwrap,'k','LineWidth',1.5);
    axis tight;
    ylabel('Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap)) max(max(OffsetRotDSUnwrap),max(meanAngBodUnwrap))], 'LineStyle','--');
    end
    colormap('parula');
    
%     for trial = 1:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(OffsetRotDS(tLimits(1),1),mean(OffsetRotDSUnwrap), ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
        cylCor = corrcoef(OffsetRotDSUnwrap,meanAngBodUnwrap);
        text(OffsetRotDS(1,1),mean(OffsetRotDSUnwrap), ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
        set(gca,'FontSize',16); 
%     end   
end

if strcmp(positionDat.exType,'LightCyl3ObjBack')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shift data for plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PVAOffsetBodAll = meanAngRaw(minFG:maxFG) - OffsetRotDS(:,2)';
    PVAOffsetObjAll = meanAngRaw(minFG:maxFG) - netAngleDS';
    for it = 1:length(meanAngRaw(minFG:maxFG))
        if PVAOffsetBodAll(it) < 0
            PVAOffsetBodAll(it) = PVAOffsetBodAll(it) + 2*pi;
        end
        if PVAOffsetObjAll(it) < 0
            PVAOffsetObjAll(it) = PVAOffsetObjAll(it) + 2*pi;
        end
    end
    PVAOffsetBodAll(isnan(PVAOffsetBodAll)) = 0;
    PVAOffsetObjAll(isnan(PVAOffsetObjAll)) = 0;
    PVAOffsetBod = mean(PVAOffsetBodAll);
    PVAOffsetObj = mean(PVAOffsetObjAll);

    % Shift the PVA by the appropriate amount 
    meanAngBod = meanAngRaw(minFG:maxFG) - PVAOffsetBod;
    meanAngBod = mod(meanAngBod-pi,2*pi)-pi;
    meanAngObj = meanAngRaw(minFG:maxFG) - PVAOffsetObj;
    meanAngObj = mod(meanAngObj-pi,2*pi)-pi;
    
    PVASpan = OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    
    % Get rid of wrap-around artifacts
    for i=1:length(OffsetRotDS)-1
         if abs(OffsetRotDS(i,2)-OffsetRotDS(i+1,2)) > pi/2
            OffsetRotDS(i,2) = NaN;
         end
         if abs(netAngle(i)-netAngle(i+1)) > pi/2
            netAngle(i) = NaN;
         end
         if abs(meanAngBod(i)-meanAngBod(i+1)) > pi/2
            meanAngBod(i) = NaN;
         end
         if abs(meanAngObj(i)-meanAngObj(i+1)) > pi/2
            meanAngObj(i) = NaN;
         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the local and global rotation to the activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    subplot(3,4,5:8);
    hold on;
    tprog = zeros(1,length(c))-pi;
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000,OffsetRotDS(:,2),'g');
    plot(PVASpan,meanAngBod,'LineWidth',1,'Color','k');
    axis tight;
    ylabel('Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(tSpan(1),0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
    
    subplot(3,4,9:12);
    hold on;
    CylWidth = 2*atan(1./sqrt(positionDat.OffsetFor.^2+positionDat.OffsetLat.^2));
    errorbarPatch2(tSpan,netAngle, CylWidth,[0 0 1]);
    plot(PVASpan,meanAngObj,'LineWidth',1,'Color','k');
    axis tight;
    ylabel('Object Orientation (rad)','FontSize',16);
    xlabel('Time (sec)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    colormap('parula');
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(netAngleDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(positionDat.tVR(1)/10000,0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
    p1=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[2*pi/3 pi pi 2*pi/3],'k');
    p2=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[-2*pi/3 -pi -pi -2*pi/3],'k');
    set(p1,'FaceAlpha',0.25); 
    set(p2,'FaceAlpha',0.25); 
    ylim([-pi pi]);
    
    netAngleH1 = OffsetRotDSUnwrap - pi/3;
    netAngleH2 = OffsetRotDSUnwrap + pi/3;
    netAngleH3 = OffsetRotDSUnwrap + 5*pi/6;
    netAngleH1 = mod(netAngleH1+pi,2*pi)-pi;
    netAngleH2 = mod(netAngleH2+pi,2*pi)-pi;
    netAngleH3 = mod(netAngleH3+pi,2*pi)-pi;

    for i=1:length(netAngleH1)-1
         if abs(netAngleH1(i)-netAngleH1(i+1)) > pi/8
            netAngleH1(i) = NaN;
         end
         if abs(netAngleH2(i)-netAngleH2(i+1)) > pi/8
            netAngleH2(i) = NaN;
         end
         if abs(netAngleH3(i)-netAngleH3(i+1)) > pi/8
            netAngleH3(i) = NaN;
         end
    end
    plot(OffsetRotDS(:,1),netAngleH1,'Color',[0 0.5 0.5],'LineWidth',2);%,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH2,'g','LineWidth',2);%,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH3,'r','LineWidth',2);%,'LineStyle',':');
    
end


if strcmp(positionDat.exType,'NoCyl3ObjBack')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shift data for plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PVAOffsetBodAll = meanAngRaw(minFG:maxFG) - OffsetRotDS(:,2)';
    PVAOffsetObjAll = meanAngRaw(minFG:maxFG) - netAngleDS';
    for it = 1:length(meanAngRaw(minFG:maxFG))
        if PVAOffsetBodAll(it) < 0
            PVAOffsetBodAll(it) = PVAOffsetBodAll(it) + 2*pi;
        end
        if PVAOffsetObjAll(it) < 0
            PVAOffsetObjAll(it) = PVAOffsetObjAll(it) + 2*pi;
        end
    end
    PVAOffsetBodAll(isnan(PVAOffsetBodAll)) = 0;
    PVAOffsetObjAll(isnan(PVAOffsetObjAll)) = 0;
    PVAOffsetBod = mean(PVAOffsetBodAll);
    PVAOffsetObj = mean(PVAOffsetObjAll);

    % Shift the PVA by the appropriate amount 
    meanAngBod = meanAngRaw(minFG:maxFG) - PVAOffsetBod;
    meanAngBod = mod(meanAngBod-pi,2*pi)-pi;
    meanAngObj = meanAngRaw(minFG:maxFG) - PVAOffsetObj;
    meanAngObj = mod(meanAngObj-pi,2*pi)-pi;
    
    PVASpan = OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    
    % Get rid of wrap-around artifacts
    for i=1:length(OffsetRotDS)-1
         if abs(OffsetRotDS(i,2)-OffsetRotDS(i+1,2)) > pi/2
            OffsetRotDS(i,2) = NaN;
         end
         if abs(netAngle(i)-netAngle(i+1)) > pi/2
            netAngle(i) = NaN;
         end
         if abs(meanAngBod(i)-meanAngBod(i+1)) > pi/2
            meanAngBod(i) = NaN;
         end
         if abs(meanAngObj(i)-meanAngObj(i+1)) > pi/2
            meanAngObj(i) = NaN;
         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the local and global rotation to the activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    subplot(3,4,5:8);
    hold on;
    tprog = zeros(1,length(c))-pi;
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000,OffsetRotDS(:,2),'g');
%     scatter(PVASpan,meanAngBod,15,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.5 0],'LineWidth',1.5);
    plot(PVASpan,meanAngBod,'Color',[0 0 0],'LineWidth',1);
    axis tight;
    ylabel('Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(tSpan(1),0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
    
    subplot(3,4,9:12);
    hold on;
%     scatter(PVASpan,meanAngObj,15,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.5 0],'LineWidth',1.5);
    plot(PVASpan,meanAngObj,'Color',[0 0 0],'LineWidth',1);
    axis tight;
    ylabel('Object Orientation (rad)','FontSize',16);
    xlabel('Time (sec)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
    end
    colormap('parula');
    
%     for trial = 2:numTrials-1
%         tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
%         cylCor = corrcoef(netAngleDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
%         text(positionDat.tVR(1)/10000,0, ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
%         set(gca,'FontSize',16); 
%     end
    
  
    p1=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[2*pi/3 pi pi 2*pi/3],'k');
    p2=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[-2*pi/3 -pi -pi -2*pi/3],'k');
    set(p1,'FaceAlpha',0.25); 
    set(p2,'FaceAlpha',0.25); 
    ylim([-pi pi]);
    
    netAngleH1 = OffsetRotDSUnwrap - pi/3;
    netAngleH2 = OffsetRotDSUnwrap + pi/3;
    netAngleH3 = OffsetRotDSUnwrap + 5*pi/6;
    netAngleH1 = mod(netAngleH1+pi,2*pi)-pi;
    netAngleH2 = mod(netAngleH2+pi,2*pi)-pi;
    netAngleH3 = mod(netAngleH3+pi,2*pi)-pi;

    for i=1:length(netAngleH1)-1
         if abs(netAngleH1(i)-netAngleH1(i+1)) > pi/8
            netAngleH1(i) = NaN;
         end
         if abs(netAngleH2(i)-netAngleH2(i+1)) > pi/8
            netAngleH2(i) = NaN;
         end
         if abs(netAngleH3(i)-netAngleH3(i+1)) > pi/8
            netAngleH3(i) = NaN;
         end
    end
    plot(OffsetRotDS(:,1),netAngleH1,'Color',[0 0.5 0.5],'LineWidth',2,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH2,'g','LineWidth',2,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH3,'r','LineWidth',2,'LineStyle',':');
    
end

%% Look at the walking S.D. vs. correlation coefficient
runCorr = zeros(length(meanAngBodUnwrap)-199,1);
runSD = zeros(length(meanAngBodUnwrap)-199,1);
for slide = 1:25:length(runCorr)
    slideCorr = corrcoef(netAngleDSUnwrap(slide:slide+199),meanAngBodUnwrap(slide:slide+199));
    runCorr(slide) = slideCorr(1,2);
    runSD(slide) = std(netAngleDSUnwrap(slide:slide+199));
end

corrSD = hist2d([runSD, runCorr], linspace(0, 3.5, 35), linspace(-1, 1, 20));

CorrVSSD = figure;
pc = pcolor(linspace(-1, 1, 19), linspace(0, 3.5, 34),corrSD);
set(pc, 'EdgeColor', 'none');
colormap(map);
axis equal;
axis tight;
xlabel('Correlation');
ylabel('Walking rotation S.D. (rad)');
set(gca,'FontSize',16); 

%% Find the correlation coefficients for the unwrapped traces

if strcmp(positionDat.exType,'LightCylFlatBack')
    unwrapPlot = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,1,1);
    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    hold on;
    tprog = zeros(1,length(c))+min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap));
    scatter(tSpan(1:10:end),tprog,3,c);
    plot(OffsetRotDS(:,1),OffsetRotDSUnwrap,'g');
    plot(OffsetRotDS(:,1),meanAngBodUnwrap,'k','LineWidth',1.5);
    axis tight;
    ylabel('Net Body Rotation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [min(min(OffsetRotDSUnwrap),min(meanAngBodUnwrap)) max(max(OffsetRotDSUnwrap),max(meanAngBodUnwrap))], 'LineStyle','--');
    end
    colormap('parula');

    for trial = 2:numTrials-1
        tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
        cylCor = corrcoef(OffsetRotDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
        text(tSpan(1),mean(OffsetRotDSUnwrap), ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
        set(gca,'FontSize',16); 

        text(tSpan(1),mean(OffsetRotDSUnwrap)/2, ['Delta: ', num2str(mean(abs(OffsetRotDSUnwrap'-mean(OffsetRotDSUnwrap)-(meanAngBodUnwrap-mean(meanAngBodUnwrap)))))],'FontSize',16);
        set(gca,'FontSize',16); 
    end
    set(gca,'FontSize',16);

    subplot(2,1,2);
    tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
    hold on;
    plot(OffsetRotDS(:,1),netAngleDSUnwrap,'b');
    plot(OffsetRotDS(:,1),meanAngBodUnwrap,'k','LineWidth',1.5);
    axis tight;
    ylabel('Net Object Orientation (rad)','FontSize',16);
    set(gca,'FontSize',16);

    for trial = 2:numTrials
        line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [min(min(netAngleDSUnwrap),min(meanAngBodUnwrap)) max(max(netAngleDSUnwrap),max(meanAngBodUnwrap))], 'LineStyle','--');
    end

    for trial = 2:numTrials-1
        tLimits = find(OffsetRotDS(:,1) >= positionDat.t(trialbreaks(trial)) & OffsetRotDS(:,1) <= positionDat.t(trialbreaks(trial+1)));
        cylCor = corrcoef(netAngleDSUnwrap(tLimits(1):tLimits(end)),meanAngBodUnwrap(tLimits(1):tLimits(end)));
        text(positionDat.tVR(1)/10000,mean(netAngleDSUnwrap), ['Corr: ', num2str(cylCor(1,2))],'FontSize',16);
        set(gca,'FontSize',16); 

        text(tSpan(1),mean(netAngleDSUnwrap)/2, ['Delta: ', num2str(mean(abs(netAngleDSUnwrap'-mean(netAngleDSUnwrap)-(meanAngBodUnwrap-mean(meanAngBodUnwrap)))))],'FontSize',16);
        set(gca,'FontSize',16); 
    end
    set(gca,'FontSize',16);
    xlabel('Time (sec)')
end

%% Find the angular difference between bump jumps and then plot when the bump jumps in the fly's trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 sec running correlation with net angle to find jumps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unwrapped, downsampled net angle = netAngleDSUnwrap
% Unwrapped activity = meanAngBodUnwrap
% Time pts = OffsetRotDS(:,1)
bumpPlotShift = input('Bump shift?');

bumpJump = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,1,2);
hold on;
% scatter(PVASpan,mod(meanAngObj+bumpPlotShift+pi,2*pi)-pi,15,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.5 0],'LineWidth',1.5);
plot(PVASpan,mod(meanAngObj+bumpPlotShift+pi,2*pi)-pi,'Color',[0 0 0],'LineWidth',1);
axis tight;
ylabel('Object Orientation (rad)','FontSize',8);
xlabel('Time (sec)','FontSize',10);
set(gca,'FontSize',10);

for trial = 2:numTrials
    line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
end
colormap('parula');

p1=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[2*pi/3 pi pi 2*pi/3],'k');
p2=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[-2*pi/3 -pi -pi -2*pi/3],'k');
set(p1,'FaceAlpha',0.25); 
set(p2,'FaceAlpha',0.25); 
ylim([-pi pi]);

netAngleH1 = OffsetRotDSUnwrap - pi/3;
netAngleH2 = OffsetRotDSUnwrap + pi/3;
netAngleH3 = OffsetRotDSUnwrap + 5*pi/6;
netAngleH1 = mod(netAngleH1+pi,2*pi)-pi;
netAngleH2 = mod(netAngleH2+pi,2*pi)-pi;
netAngleH3 = mod(netAngleH3+pi,2*pi)-pi;

for i=1:length(netAngleH1)-1
     if abs(netAngleH1(i)-netAngleH1(i+1)) > pi/8
        netAngleH1(i) = NaN;
     end
     if abs(netAngleH2(i)-netAngleH2(i+1)) > pi/8
        netAngleH2(i) = NaN;
     end
     if abs(netAngleH3(i)-netAngleH3(i+1)) > pi/8
        netAngleH3(i) = NaN;
     end
end
    plot(OffsetRotDS(:,1),netAngleH1,'Color',[0 0.5 0.5],'LineWidth',2,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH2,'g','LineWidth',2,'LineStyle',':');
    plot(OffsetRotDS(:,1),netAngleH3,'r','LineWidth',2,'LineStyle',':');
if strcmp(positionDat.exType,'LightCyl3ObjBack')
    errorbarPatch2(tSpan,netAngle, CylWidth,[0 0 1]);
end

% Find the running correlation
tStep = mean(diff(OffsetRotDS(:,1)));
tCorr = 5; % # of secs for running corr.
nPts = round(tCorr/tStep);

runCorr = zeros(length(OffsetRotDS)-nPts+1,1);
for runIt = 1:length(OffsetRotDS)-nPts+1
    matCorr = corrcoef(netAngleDSUnwrap(runIt:runIt+nPts-1),meanAngBodUnwrap(runIt:runIt+nPts-1));
    runCorr(runIt) = matCorr(1,2);
end

subplot(4,1,3);
hold on;
plot(OffsetRotDS(1:end-nPts+1,1)+round(tCorr/2),runCorr);
set(gca,'FontSize',10);
xlabel('Time (sec)');
ylabel(['Running Correlation (',num2str(tCorr), ' sec)'],'FontSize',8);
axis tight;
xlim([OffsetRotDS(1,1) OffsetRotDS(end,1)])
line([OffsetRotDS(1,1) OffsetRotDS(end,1)], [0 0],'Color',[0 0 0]);

for trial = 2:numTrials
    line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-1 1], 'LineStyle','--');
end
xlim([positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000]);

% [jump,jumpLocs] = findpeaks(-runCorr,OffsetRotDS(1:end-nPts+1,1),'MinPeakDistance',tCorr);
% jumpLoc = jumpLocs(find(jump>-0.25));
% jump = -jump(find(jump>-0.25));
% scatter(jumpLoc,jump);
% [locked,lockedLocs] = findpeaks(runCorr,OffsetRotDS(1:end-nPts+1,1),'MinPeakDistance',tCorr);
% scatter(lockedLocs(find(locked>0)),locked(find(locked>0)));

% Find the running difference
tDiff = 5; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiff = zeros(length(OffsetRotDS)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotDS)-nPtsDiff+1
    runDiff(runIt) = mean(abs(netAngleDSUnwrap(runIt:runIt+nPtsDiff-1)'-netAngleDSUnwrap(runIt)'-(meanAngBodUnwrap(runIt:runIt+nPtsDiff-1)-meanAngBodUnwrap(runIt))));
end

subplot(4,1,4);
hold on;
plot(OffsetRotDS(1:end-nPtsDiff+1,1)+round(tDiff/2),runDiff);
set(gca,'FontSize',10);
xlabel('Time (sec)');
ylabel(['Running Difference (',num2str(tDiff), ' sec)'],'FontSize',8);
axis tight;
xlim([OffsetRotDS(1,1) OffsetRotDS(end,1)])
[jumpDiff,jumpDiffLocs] = findpeaks(runDiff,OffsetRotDS(1:end-nPtsDiff+1,1),'MinPeakDistance',5);
jumpDiffLocs = jumpDiffLocs(find(jumpDiff>3*pi/8));
jumpDiff = jumpDiff(find(jumpDiff>3*pi/8));
scatter(jumpDiffLocs+round(tDiff/2),jumpDiff);
badPks = input('List bad peaks[#]');
jumpDiffLocs(badPks) = [];
jumpDiff(badPks) = [];
xlim([positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000]);

jumpAmts = zeros(length(jumpDiff),1);

for jumps = 1:length(jumpDiff)
    subplot(4,1,2);
    line([jumpDiffLocs(jumps)+tDiff/2 jumpDiffLocs(jumps)+tDiff/2],[ -pi pi],'Color',[0 0 0],'LineStyle','--');
    if jumps > 1 & jumps < length(jumpDiff)
        prevMax = find(runDiff == ...
            min(runDiff(...
            find(OffsetRotDS == jumpDiffLocs(jumps-1)):...
            find(OffsetRotDS == jumpDiffLocs(jumps)))))...
            +round(nPtsDiff/2);
        nextMax = find(runDiff...
            == min(runDiff(...
            find(OffsetRotDS == jumpDiffLocs(jumps)):...
            find(OffsetRotDS == jumpDiffLocs(jumps+1)))))...
            +round(nPtsDiff/2);
    elseif jumps == 1
        prevMax = find(runDiff...
            == min(runDiff(...
            1:find(OffsetRotDS == jumpDiffLocs(jumps)))))...
            +round(nPtsDiff/2);
        nextMax = find(runDiff...
            == min(runDiff(....
            find(OffsetRotDS == jumpDiffLocs(jumps))...
            :find(OffsetRotDS == jumpDiffLocs(jumps+1)))))...
            +round(nPtsDiff/2);
    else
        prevMax = find(runDiff...
            == min(runDiff(...
            find(OffsetRotDS == jumpDiffLocs(jumps-1)):...
            find(OffsetRotDS == jumpDiffLocs(jumps)))))...
            +round(nPtsDiff/2);
        nextMax = find(runDiff...
            == min(runDiff(...
            find(OffsetRotDS == jumpDiffLocs(jumps)):end)));
    end
%     rectangle('Position',[OffsetRotDS(prevMax), -pi, OffsetRotDS(nextMax) - OffsetRotDS(prevMax)+tDiff, 2*pi],'EdgeColor','b');
    text(jumpDiffLocs(jumps)+tDiff/2, 11*pi/9,num2str(jumps),'FontSize',16);
    xlim([positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000]);
    
    subplot(4,1,1)
    slop = runDiff(nextMax-round(nPtsDiff/2)) - runDiff(prevMax-round(nPtsDiff/2));
    angMov = mean(OffsetRotDSUnwrap(nextMax:nextMax+nPtsDiff-1)) - mean(OffsetRotDSUnwrap(prevMax:prevMax+nPtsDiff-1));
    bumpMov = mean(meanAngBodUnwrap(nextMax:nextMax+nPtsDiff-1))...
        - mean(meanAngBodUnwrap(prevMax:prevMax+nPtsDiff-1));
    if jumps ==1
        text(jumps-1, 4,'Jump ','FontSize',16);
        text(jumps-1, 3,'Slop: ','FontSize',8);
        text(jumps-1, 2,'Body shift: ','FontSize',8);
        text(jumps-1, 1,'Bump shift: ','FontSize',8);
        text(jumps-1, 0,'Net shift: ','FontSize',8);
    end 
    text(jumps, 4,num2str(jumps),'FontSize',16);
    text(jumps, 3,[num2str(slop/pi,3),'*\pi'],'FontSize',8);
    text(jumps, 2,[num2str(angMov/pi,3),'*\pi'],'FontSize',8);
    text(jumps, 1,[num2str(bumpMov/pi,3),'*\pi'],'FontSize',8);
    text(jumps, 0,[num2str((angMov - bumpMov)/pi,3),'*\pi'],'FontSize',8);
    xlim([1 length(jumpDiff)+1]);
    ylim([0 4]);
    axis off;
    
    jumpAmts(jumps) = (angMov - bumpMov)/pi;
    
end


histJumps = figure;
edges = [-1/24:1/12:2+1/24];
jumpHist = histogram(jumpAmts,edges);

xlim([0 2]);
line([0.5 0.5], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
line([2/3 2/3], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
line([5/6 5/6], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
line([0.5+2/3 0.5+2/3], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
line([0.5+5/6 0.5+5/6], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
line([2/3+5/6 2/3+5/6], [0 0.5],'LineStyle','--','Color',[0 0 0],'LineWidth',1.5);
set(gca,'FontSize',16);
xlabel('Jump Amount (*\pi rad)');
ylabel('Counts');

bumpJumpPos = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
c = colormap(hsv(length(positionDat.OffsetLat)));
scatter(positionDat.OffsetLat,positionDat.OffsetFor,1,c);
for jumpB = 1:length(jumpDiffLocs)
    jumpPt = find(positionDat.t == jumpDiffLocs(jumpB));
    scatter(positionDat.OffsetLat(jumpPt),positionDat.OffsetFor(jumpPt),60,[0 0 0],'o','LineWidth',3);
end

 if (strcmp(positionDat.exType,'LightCylFlatBack') || strcmp(positionDat.exType,'LightCyl3ObjBack'))
    viscircles([0 0], 1,'EdgeColor','b');
end
if (strcmp(positionDat.exType,'NoCyl3ObjBack') || strcmp(positionDat.exType,'LightCyl3ObjBack'))
    line([0 2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 0 0]);
    line([0 -2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 1 0]);
    line([0 2*cos(pi/6-pi/2)],[0 2*sin(pi/6-pi/2)],'LineStyle','--','Color',[1 0 0]);
end
axis equal;
xlabel('Distance (cm)');
ylabel('Distance (cm)');
set(gca,'FontSize',16);
scatter(0,-5,100,'p','k');

%% Stick the bump to the local object and see what happens
% Find the running correlation
tStep = mean(diff(OffsetRotDS(:,1)));
tCorr = 5; % # of secs for running corr.
nPts = round(tCorr/tStep);

runCorrGlob = zeros(length(OffsetRotDS)-nPts+1,1);
runCorrObj = zeros(length(OffsetRotDS)-nPts+1,1);
for runIt = 1:length(OffsetRotDS)-nPts+1
    matCorrGlob = corrcoef(OffsetRotDSUnwrap(runIt:runIt+nPts-1),meanAngBodUnwrap(runIt:runIt+nPts-1));
    runCorrGlob(runIt) = matCorrGlob(1,2);
    matCorrObj = corrcoef(netAngleDSUnwrap(runIt:runIt+nPts-1),meanAngBodUnwrap(runIt:runIt+nPts-1));
    runCorrObj(runIt) = matCorrObj(1,2);
end


% Find the running difference
tDiff = 10; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiffGlob = zeros(length(OffsetRotDS)-nPtsDiff+1,1);
runDiffObj = zeros(length(OffsetRotDS)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotDS)-nPtsDiff+1
    runDiffGlob(runIt) = mean(abs(OffsetRotDSUnwrap(runIt:runIt+nPtsDiff-1)'-OffsetRotDSUnwrap(runIt)-(meanAngBodUnwrap(runIt:runIt+nPtsDiff-1)-meanAngBodUnwrap(runIt))));
    runDiffObj(runIt) = mean(abs(netAngleDSUnwrap(runIt:runIt+nPtsDiff-1)'-netAngleDSUnwrap(runIt)-(meanAngBodUnwrap(runIt:runIt+nPtsDiff-1)-meanAngBodUnwrap(runIt))));
end


% Plot the metrics
locVSglob = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,1,1)
hold on;
plot(OffsetRotDS(1:end-nPts+1,1),runCorrGlob,'g');
plot(OffsetRotDS(1:end-nPts+1,1),runCorrObj,'b');
set(gca,'FontSize',14);
title(['Running Correlation (',num2str(tCorr), ' sec)']);
axis tight;
xlim([OffsetRotDS(1,1) OffsetRotDS(end,1)])

subplot(3,1,2)
hold on;
plot(OffsetRotDS(1:end-nPtsDiff+1,1),runDiffGlob,'g');
plot(OffsetRotDS(1:end-nPtsDiff+1,1),runDiffObj,'b');
set(gca,'FontSize',14);
title(['Running Difference (',num2str(tDiff), ' sec)']);
ylabel('Radians');
axis tight;
xlim([OffsetRotDS(1,1) OffsetRotDS(end,1)])


% Find where the PVA most closely follows the local object 
stuckLocal = find((runDiffGlob-runDiffObj) == max(runDiffGlob-runDiffObj));
line([OffsetRotDS(stuckLocal,1) OffsetRotDS(stuckLocal,1)],[0 max(runDiffGlob)],'Color',[0 0 0],'LineStyle','--');

subplot(3,1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift data for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PVASpan = OffsetRotDS(:,1)-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;

OffsetRotDSUnwrapPlot = mod(OffsetRotDSUnwrap+pi,2*pi)-pi;
meanAngBodUnwrapPlot = mod(meanAngBodUnwrap-meanAngBodUnwrap(stuckLocal+round(nPtsDiff/2))+netAngleDSUnwrap(stuckLocal+round(nPtsDiff/2))+pi,2*pi)-pi;
netAngleDSUnwrapPlot = mod(netAngleDSUnwrap+pi,2*pi)-pi;

% Get rid of wrap-around artifacts
for i=1:length(OffsetRotDS)-1
     if abs(netAngleDSUnwrapPlot(i)-netAngleDSUnwrapPlot(i+1)) > pi/2
        netAngleDSUnwrapPlot(i) = NaN;
     end
     if abs(meanAngBodUnwrapPlot(i)-meanAngBodUnwrapPlot(i+1)) > pi/2
        meanAngBodUnwrapPlot(i) = NaN;
     end
end

hold on;
tSpan = positionDat.t-positionDat.t(1)+(positionDat.tVR(1)-positionDat.tFrameGrab(1))/10000;
CylWidth = 2*atan(1./sqrt(positionDat.OffsetFor.^2+positionDat.OffsetLat.^2));
errorbarPatch2(tSpan,netAngle, CylWidth,[0 0 1]);
plot(PVASpan,meanAngBodUnwrapPlot,'LineWidth',1,'Color','k');
axis tight;
set(gca,'FontSize',14);
title('Object Orientation vs. PVA ');
ylabel('Radians')
xlabel('Time (sec)');

for trial = 2:numTrials
    line([positionDat.t(trialbreaks(trial)) positionDat.t(trialbreaks(trial))], [-pi pi], 'LineStyle','--');
end

p1=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[2*pi/3 pi pi 2*pi/3],'k');
p2=patch([positionDat.tVR(1)/10000 positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000 positionDat.t(end)-positionDat.t(1)+positionDat.tVR(1)/10000],[-2*pi/3 -pi -pi -2*pi/3],'k');
set(p1,'FaceAlpha',0.25); 
set(p2,'FaceAlpha',0.25); 
ylim([-pi pi]);

netAngleH1 = OffsetRotDSUnwrap - pi/3;
netAngleH2 = OffsetRotDSUnwrap + pi/3;
netAngleH3 = OffsetRotDSUnwrap + 5*pi/6;
netAngleH1 = mod(netAngleH1+pi,2*pi)-pi;
netAngleH2 = mod(netAngleH2+pi,2*pi)-pi;
netAngleH3 = mod(netAngleH3+pi,2*pi)-pi;

for i=1:length(netAngleH1)-1
     if abs(netAngleH1(i)-netAngleH1(i+1)) > pi/8
        netAngleH1(i) = NaN;
     end
     if abs(netAngleH2(i)-netAngleH2(i+1)) > pi/8
        netAngleH2(i) = NaN;
     end
     if abs(netAngleH3(i)-netAngleH3(i+1)) > pi/8
        netAngleH3(i) = NaN;
     end
end
plot(OffsetRotDS(:,1),netAngleH1,'Color',[0 0.5 0.5],'LineWidth',1);%,'LineStyle',':');
plot(OffsetRotDS(:,1),netAngleH2,'g','LineWidth',1);%,'LineStyle',':');
plot(OffsetRotDS(:,1),netAngleH3,'r','LineWidth',1);%,'LineStyle',':');

% Look at how much time the bump sticks to the local object

% Downsample the cylinder width
CylWidth = 2*atan(1./sqrt(positionDat.OffsetFor.^2+positionDat.OffsetLat.^2));
minFG = ceil(max(find(positionDat.tFrameGrab <= positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);
CylWidthDS = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tDS = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    CylWidthDS(interp-minFG+1) = CylWidth(tDS(1));
end

% Find the difference between the object centers and the PVA
diffObj = zeros(length(OffsetRotDS),1);
diffH1 = zeros(length(OffsetRotDS),1);
diffH2 = zeros(length(OffsetRotDS),1);
diffH3 = zeros(length(OffsetRotDS),1);

diffObj = min(abs(netAngleDS-meanAngBodUnwrapPlot'),2*pi-abs(netAngleDS-meanAngBodUnwrapPlot'));
diffH1 = min(abs(netAngleH1-meanAngBodUnwrapPlot'),2*pi-abs(netAngleH1-meanAngBodUnwrapPlot'));
diffH2 = min(abs(netAngleH2-meanAngBodUnwrapPlot'),2*pi-abs(netAngleH2-meanAngBodUnwrapPlot'));
diffH3 = min(abs(netAngleH3-meanAngBodUnwrapPlot'),2*pi-abs(netAngleH3-meanAngBodUnwrapPlot'));

% Figure out which feature the bump is locked to.
% Object == 0, H1 == 1, H2 == 2, H3 == 3
bumpLockID = zeros(length(OffsetRotDS),1);
for IDassign = 1:length(bumpLockID)
    if diffObj(IDassign) < CylWidthDS(IDassign)
        continue;
    elseif (diffObj(IDassign)-CylWidth(IDassign)) < min([diffH1(IDassign),diffH2(IDassign),diffH3(IDassign)])
        continue;
    elseif diffH1(IDassign) < min([diffH2(IDassign),diffH3(IDassign)])
        bumpLockID(IDassign) = 1;
    elseif diffH2(IDassign) < min([diffH1(IDassign),diffH3(IDassign)])
        bumpLockID(IDassign) = 2;
    elseif diffH3(IDassign) < min([diffH1(IDassign),diffH2(IDassign)])
        bumpLockID(IDassign) = 3;
    end
end
figure;
hold on;
scatter(PVASpan(find(bumpLockID == 0)),netAngleDS(find(bumpLockID == 0)) + CylWidthDS(find(bumpLockID == 0)),10,'filled','b');
scatter(PVASpan(find(bumpLockID == 0)),netAngleDS(find(bumpLockID == 0)) - CylWidthDS(find(bumpLockID == 0)),10,'filled','b');
scatter(PVASpan(find(bumpLockID == 1)),netAngleH1(find(bumpLockID == 1)),10,[0 0.5 0.5],'filled');
scatter(PVASpan(find(bumpLockID == 2)),netAngleH2(find(bumpLockID == 2)),10,'filled','g');
scatter(PVASpan(find(bumpLockID == 3)),netAngleH3(find(bumpLockID == 3)),10,'filled','r');
plot(PVASpan,meanAngBodUnwrapPlot,'LineWidth',1,'Color','k');
axis tight;
set(gca,'FontSize',14);
ylabel('Radians')
xlabel('Time (sec)');


length(find(bumpLockID == 0))/length(bumpLockID)

%% Look at the periods where the local object is invisible

% Downsample the cylinder width
CylWidth = 2*atan(1./sqrt(positionDat.OffsetFor.^2+positionDat.OffsetLat.^2));
minFG = ceil(max(find(positionDat.tFrameGrab <= positionDat.tVR(1)))/num_planes);
maxFG = round(length(positionDat.tFrameGrab)/num_planes);
CylWidthDS = zeros(maxFG-minFG+1,1);
for interp = minFG:maxFG
    tDS = find(positionDat.t >= (positionDat.t(1) + (positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
    CylWidthDS(interp-minFG+1) = CylWidth(tDS(1));
end

% Find the regions where the cylinder is not visible

% Find the running difference
tStep = mean(diff(OffsetRotDS(:,1)));
tDiff = 10; % # of secs for running corr.
nPtsDiff = round(tDiff/tStep);

runDiff = zeros(length(OffsetRotDS)-nPtsDiff+1,1);
for runIt = 1:length(OffsetRotDS)-nPtsDiff+1
    runDiff(runIt) = mean(abs(netAngleDSUnwrap(runIt:runIt+nPtsDiff-1)'-netAngleDSUnwrap(runIt)-(meanAngBodUnwrap(runIt:runIt+nPtsDiff-1)-meanAngBodUnwrap(runIt))));
end


netAngle2 = mod(netAngleDSUnwrap,2*pi)-pi;
invis = find(((netAngle2+CylWidthDS) < pi*1/3) & ((netAngle2-CylWidthDS) > -pi*1/3));
invis = invis(find(invis<length(runDiff)+round(nPtsDiff/2)));

% Find the region bounds
invisBreaks = find(abs(diff(invis))>=2);
invisBreaks = vertcat(invisBreaks, invisBreaks+1);
invisBreaks = vertcat(invisBreaks, 1, length(invis));
invisBreaks = sortrows(invisBreaks);

% Find where the running difference corresponds to the invisible object and
% is less than pi/4
bumpLock = find(runDiff(invis-round(nPtsDiff/2))<pi/4);

% Plot the regions
topEdge = mod(netAngleDSUnwrap+CylWidthDS-pi,2*pi)-pi;
bottomEdge = mod(netAngleDSUnwrap-CylWidthDS-pi,2*pi)-pi;

% Plot interesting things
objInvis = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
hold on;
plot(PVASpan,meanAngObj,'k');
scatter(OffsetRotDS(:,1),topEdge,4,'b','filled');
scatter(OffsetRotDS(:,1),bottomEdge,4,'b','filled');
scatter(OffsetRotDS(invis,1),topEdge(invis),8,'r');
scatter(OffsetRotDS(invis,1),bottomEdge(invis),8,'r');
% for brk = 1:length(invisBreaks)
%     line([OffsetRotDS(invis(invisBreaks(brk)),1) OffsetRotDS(invis(invisBreaks(brk)),1)],[-pi pi]);
% end
xlim([OffsetRotDS(1,1) OffsetRotDS(end-nPtsDiff+1,1)]);
ylim([-pi pi]);
line([OffsetRotDS(1,1) OffsetRotDS(end-nPtsDiff+1,1)],[2*pi/3 2*pi/3],'Color','k','LineStyle','--');
line([OffsetRotDS(1,1) OffsetRotDS(end-nPtsDiff+1,1)],[-2*pi/3 -2*pi/3],'Color','k','LineStyle','--');
set(gca,'FontSize',16)
ylabel('Position (rad)');

subplot(2,1,2);
hold on;
plot(OffsetRotDS(1:end-nPtsDiff+1,1)+OffsetRotDS(round(nPtsDiff/2),1),runDiff);
scatter(OffsetRotDS(invis-round(nPtsDiff/2),1)+OffsetRotDS(round(nPtsDiff/2),1),runDiff(invis-round(nPtsDiff/2)),'r');
scatter(OffsetRotDS(invis(bumpLock)-round(nPtsDiff/2),1)+OffsetRotDS(round(nPtsDiff/2),1),runDiff(invis(bumpLock)-round(nPtsDiff/2)),'g');
xlim([OffsetRotDS(1,1) OffsetRotDS(end-nPtsDiff+1,1)]);
line([OffsetRotDS(1,1) OffsetRotDS(end-nPtsDiff+1,1)],[pi/4 pi/4],'Color','k','LineStyle','--');
set(gca,'FontSize',16)
ylabel('Running difference (rad)');
xlabel('Time (sec)');

text(OffsetRotDS(end-nPtsDiff+1,1)/2, max(runDiff)*9/10, ['% unseen: ', num2str(length(invis)/length(OffsetRotDS)*100)],'FontSize',24);
text(OffsetRotDS(end-nPtsDiff+1,1)/2, max(runDiff)*8/10, ['% unseen locked: ', num2str(length(bumpLock)/length(invis)*100)],'FontSize',24);

%% Save the data
set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(act,strcat(fullpath{1}(1:end-4),'_BumpAct'),'-dpdf'); 
set(rotComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(rotComp,strcat(fullpath{1}(1:end-4),'_BumpRot'),'-dpdf'); 
if strcmp(positionDat.exType,'NoCyl3ObjBack') || strcmp(positionDat.exType,'NoCylFlatBack')
    set(CorrVSSD,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(CorrVSSD,strcat(fullpath{1}(1:end-4),'_CorrSD'),'-dpdf'); 
end
if strcmp(positionDat.exType,'LightCylFlatBack')
    set(unwrapPlot,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(unwrapPlot,strcat(fullpath{1}(1:end-4),'_Unwrap'),'-dpdf'); 
    set(objInvis,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(objInvis,strcat(fullpath{1}(1:end-4),'_objInvis'),'-dpdf'); 
end
if strcmp(positionDat.exType,'NoCyl3ObjBack')
    set(bumpJump,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(bumpJump,strcat(fullpath{1}(1:end-4),'_BumpJump'),'-dpdf'); 
    set(bumpJumpPos,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(bumpJumpPos,strcat(fullpath{1}(1:end-4),'_BumpJumpLoc'),'-dpdf');
    set(histJumps,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(histJumps,strcat(fullpath{1}(1:end-4),'_BumpJumpHist'),'-dpdf');
end
if strcmp(positionDat.exType,'LightCyl3ObjBack')
    set(locVSglob,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(locVSglob,strcat(fullpath{1}(1:end-4),'_BumpJump'),'-dpdf'); 
end

save(strcat(fullpath{1}(1:end-10),'.mat'),'fullpath','ROIs','ROIaveREFMax','ROIaveREFMean','positionDat');
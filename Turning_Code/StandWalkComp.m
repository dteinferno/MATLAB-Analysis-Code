%% Get out the walking and standing data
% Set parameters for extracting walking bouts
tMergeStand = 0.01;
idxFilt = 1;
tMerge = 0.8;
tLength = 1.5;

standingDat = {};
walkingDat = {};

for flyID = 1:numFlies
    for darkID = 1:1:length(allFlyData{flyID}.dark)
        if flyID == 1 & (darkID == 6 || darkID == 7)
            continue;
        end
        if flyID == 4 & (darkID == 1 || darkID == 2 || darkID == 3)
            continue;
        end
        if flyID == 10 & (darkID == 5)
            continue;
        end

        % Match the behavior to the imaging
        num_planes = round(length(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab)/...
            length(allFlyData{flyID}.dark{darkID}.RROIaveMax));
        tSpan = allFlyData{flyID}.dark{darkID}.positionDat.t - ...
            allFlyData{flyID}.dark{darkID}.positionDat.t(1) + ...
            allFlyData{flyID}.dark{darkID}.positionDat.tVR(1)/10000;
        minFG = ceil(min(find(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab >=...
            allFlyData{flyID}.dark{darkID}.positionDat.tVR(1)))/num_planes);
        maxFG = round(length(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab)/num_planes);

        % Match the position data to the framegrab times
        OffsetRotMatch = zeros(maxFG-minFG+1,2);
        OffsetForMatch = zeros(maxFG-minFG+1,1);
        OffsetLatMatch = zeros(maxFG-minFG+1,1);
        for interp = minFG:maxFG
            tMatch = find(allFlyData{flyID}.dark{darkID}.positionDat.t >=...
                (allFlyData{flyID}.dark{darkID}.positionDat.t(1) +...
                (allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                allFlyData{flyID}.dark{darkID}.positionDat.tVR(1))/10000));
            OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.dark{darkID}.positionDat.t(tMatch(1));
            OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.dark{darkID}.positionDat.OffsetRot(tMatch(1));
            OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.dark{darkID}.positionDat.OffsetFor(tMatch(1));
            OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.dark{darkID}.positionDat.OffsetLat(tMatch(1));
        end

        OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);

        % Find walking and standing bouts
        clear motiondata;
        motiondata{1}(:,1) = OffsetRotMatch(1:end-1,1);
        motiondata{1}(:,2) = diff(OffsetRotMatchUnwrap);
        motiondata{1}(:,3) = sqrt(diff(OffsetForMatch).^2+diff(OffsetLatMatch).^2);

        [upWalkIdxs, downWalkIdxs, walkIdxs, startStandIdxs, endStandIdxs] = extractWalkingBoutInfo(motiondata, idxFilt, tMerge, tMergeStand, tLength);

        % Savitzky-Golay Filter the activity
        sgolayOrder = 3;
        sgolayWindow = 7;
        GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMax(:,minFG:maxFG),sgolayOrder,sgolayWindow,[],2);

        % Find PVA of the activity
        num_ROIs = 16;
        angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
        angsraw = angsraw';
        clear meanAngRaw;
        for ts = 1:length(GROIaveMaxFilt)
            meanAngRaw(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
        end

        % figure;
        % subplot(2,1,1);
        % imagesc(GROIaveMaxFilt);
        % hold on;
        % plot(8/pi*(meanAngRaw+pi));

        % For the standing bouts, remove the bump ROIs
        for stnd = 1:length(startStandIdxs{1})
            PVApk = round(8/pi*(meanAngRaw(startStandIdxs{1}(stnd)) + pi));
            PVArng = [PVApk-2:PVApk+2];
            PVArng(find(PVArng<1)) = PVArng(find(PVArng<1))+16;
            PVArng(find(PVArng>16)) = PVArng(find(PVArng>16))-16;
            GROIaveMaxFilt(PVArng,startStandIdxs{1}(stnd):endStandIdxs{1}(stnd)) = zeros(length(PVArng),endStandIdxs{1}(stnd)-startStandIdxs{1}(stnd)+1);;
        end

        % For the walking bouts, remove the bump ROIs
        for wlk = 1:length(upWalkIdxs{1})
            for step =1:downWalkIdxs{1}(wlk)-upWalkIdxs{1}(wlk)+1
                PVApk = round(8/pi*(meanAngRaw(upWalkIdxs{1}(wlk)+step-1) + pi));
                PVArng = [PVApk-2:PVApk+2];
                PVArng(find(PVArng<1)) = PVArng(find(PVArng<1))+16;
                PVArng(find(PVArng>16)) = PVArng(find(PVArng>16))-16;
                GROIaveMaxFilt(PVArng,upWalkIdxs{1}(wlk)+step-1) = zeros(length(PVArng),1);
            end
        end

        % subplot(2,1,2);
        % imagesc(GROIaveMaxFilt);
        % hold on;
        % plot(8/pi*(meanAngRaw+pi));

        % Collapse the ROIs onto one another
        GROIaveMaxFiltCollapse = mean(GROIaveMaxFilt,1)*16/11;

        for stnd = 1:length(startStandIdxs{1})
            standingDat{flyID}{darkID}{stnd} = GROIaveMaxFiltCollapse(startStandIdxs{1}(stnd):endStandIdxs{1}(stnd));
        end
        for wlk = 1:length(upWalkIdxs{1})
            walkingDat{flyID}{darkID}{wlk} = GROIaveMaxFiltCollapse(upWalkIdxs{1}(wlk):downWalkIdxs{1}(wlk));
        end
    end
end

%% Collapse all of the various data points together
framerate = 1./mean(diff(OffsetRotMatch(:,1)));

tStartStand = round( 3 * framerate);
tStopStand = round( 3 * framerate);
tStartWalk = round( 3 * framerate);
tStopWalk = round( 3 * framerate);

standBegin = zeros(numFlies,tStartStand);
standEnd = zeros(numFlies,tStopStand);
walkBegin = zeros(numFlies,tStartWalk);
walkEnd = zeros(numFlies,tStopWalk);

totalStandBegin = 0;
totalStandEnd = 0;
totalWalkBegin = 0;
totalWalkEnd = 0;

for stnd=1:length(standingDat)
    for trial = 1:length(standingDat{stnd})
        for bout = 1:length(standingDat{stnd}{trial})
            if length(standingDat{stnd}{trial}{bout}) >= tStartStand
                standBegin(stnd,:) = standBegin(stnd,:) + standingDat{stnd}{trial}{bout}(1:tStartStand);
                totalStandBegin = totalStandBegin + 1;
            end
            if length(standingDat{stnd}{trial}{bout}) >= tStopStand
                standEnd(stnd,:) = standEnd(stnd,:) + standingDat{stnd}{trial}{bout}(end-tStopStand+1:end);
                totalStandEnd = totalStandEnd + 1;
            end
        end 
    end
    standBegin(stnd,:) = standBegin(stnd,:)./totalStandBegin;
    totalStandBegin = 0;
    standEnd(stnd,:) = standEnd(stnd,:)./totalStandEnd;
    totalStandEnd = 0;
end

for wlk=1:length(walkingDat)
    for trial = 1:length(walkingDat{wlk})
        for bout = 1:length(walkingDat{wlk}{trial})
            if length(walkingDat{wlk}{trial}{bout}) >= tStartWalk
                walkBegin(wlk,:) = walkBegin(wlk,:) + walkingDat{wlk}{trial}{bout}(1:tStartWalk);
                totalWalkBegin = totalWalkBegin + 1;
            end
            if length(walkingDat{wlk}{trial}{bout}) >= tStopWalk
                walkEnd(wlk,:) = walkEnd(wlk,:) + walkingDat{wlk}{trial}{bout}(end-tStopWalk+1:end);
                totalWalkEnd = totalWalkEnd + 1;
            end
        end
    end
    walkBegin(wlk,:) = walkBegin(wlk,:)./totalWalkBegin;
    totalWalkBegin = 0;
    walkEnd(wlk,:) = walkEnd(wlk,:)./totalWalkEnd;
    totalWalkEnd = 0;
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,2,1);
hold on;
plot(OffsetRotMatch(1:length(standBegin))-OffsetRotMatch(1),standBegin');
plot(OffsetRotMatch(1:length(standBegin))-OffsetRotMatch(1),mean(standBegin,1),'Color','k','LineWidth',2);
set(gca,'FontSize',16);
title('Beginning of standing bout');
ylim([1 2]);

subplot(4,2,3);
hold on;
plot(OffsetRotMatch(1:length(standEnd))-OffsetRotMatch(1),standEnd');
plot(OffsetRotMatch(1:length(standEnd))-OffsetRotMatch(1),mean(standEnd,1),'Color','k','LineWidth',2);
set(gca,'FontSize',16);
title('End of standing bout');
ylim([1 2]);

subplot(4,2,5);
hold on;
plot(OffsetRotMatch(1:length(walkBegin))-OffsetRotMatch(1),walkBegin');
plot(OffsetRotMatch(1:length(walkBegin))-OffsetRotMatch(1),mean(walkBegin,1),'Color','k','LineWidth',2);
set(gca,'FontSize',16);
title('Beginning of walking bout');
ylim([1 2]);

subplot(4,2,7);
hold on;
plot(OffsetRotMatch(1:length(walkEnd))-OffsetRotMatch(1),walkEnd');
plot(OffsetRotMatch(1:length(walkEnd))-OffsetRotMatch(1),mean(walkEnd,1),'Color','k','LineWidth',2);
set(gca,'FontSize',16);
title('End of walking bout');
xlabel('Time (sec)');
ylim([1 2]);
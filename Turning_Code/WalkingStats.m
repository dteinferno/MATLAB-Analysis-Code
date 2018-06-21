%% Load all of the EB data
% Specify directories
dirRoot = 'D:\Imaging\2Color';
region = 'EB';

imDir{1} = '20160504';
redLine = '37F06';
greenLine = '60D05';

% 
% imDir{1} = '20160423-Flip';
% imDir{2} = '20160425-Flip';
% imDir{3} = '20160426-Flip';
% imDir{4} = '20160429';
% imDir{5} = '20160505';
% redLine = '60D05';
% greenLine = '37F06';

% Get the fly IDs
allFlyData{1}.ID = 'Empty';
numFlies = 1;
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            for flyStep = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyStep}.ID)
                    break;
                end
                if flyStep == length(allFlyData)
                    allFlyData{numFlies}.ID = flyName;
                    numFlies = numFlies + 1;
                end
            end
        end
    end
end
numFlies = numFlies-1;

% Get the position data and trial ID for all of the flies and all of the trials
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs') & ~strcmpi(fileName(end-5:end-4),'NO') & ~strcmpi(fileName(end-8:end-4),'NORaw') 
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
            
            % Do some more preprocessing
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 11;
            RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(GROIaveMean,sgolayOrder,sgolayWindow,[],2);
            
            % Match the behavior to the imaging
            num_planes = round(length(positionDat.tFrameGrab)/length(RROIaveMax));
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
            vF = sqrt(diff(OffsetForMatch).^2+diff(OffsetLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));
            vF = sgolayfilt(vF,sgolayOrder,sgolayWindow);
            vSS = diff(OffsetLatMatch)./mean(diff(OffsetRotMatch(:,1)));
            
            positionDatMatch = {};
            positionDatMatch.OffsetRotMatch = OffsetRotMatch;
            positionDatMatch.OffsetForMatch = OffsetForMatch;
            positionDatMatch.OffsetLetMatch = OffsetLatMatch;
            positionDatMatch.vRot = vRot;
            positionDatMatch.vF = vF;
            positionDatMatch.vSS = vSS;
            positionDatMatch.minFG = minFG;
            positionDatMatch.maxFG = maxFG;
            
            for flyID = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyID}.ID)
                    fileParts = strsplit(fileName,'_');
                         if strcmpi(fileParts(5),'1x')
                            allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'2x')
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'Dark')
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);  
                    end
                end
            end
        end 
    end
end

%% Load all of the PB data
% Specify directories
dirRoot = 'D:\Imaging\2Color';
region = 'PB';

% imDir{1} = '20160621';
% imDir{2} = '20160908';
% redLine = '37F06';
% redSpan = [1:8 11:18];
% greenLine = '60D05';
% greenSpan = [2:17];

% imDir{1} = '20160622';
% imDir{2} = '20160907';
% imDir{3} = '20160908-Flip';
% redLine = '60D05';
% redSpan = [2:17];
% greenLine = '37F06';
% greenSpan = [1:8 11:18];

imDir{1} = '20160903';
redLine = 'VT8135';
redSpan = [1:8 11:18];
greenLine = '60D05';
greenSpan = [2:17];

% Get the fly IDs
allFlyData{1}.ID = 'Empty';
numFlies = 1;
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs') & ~strcmpi(fileName(1:3),'PVA') & ~strcmpi(fileName(1:3),'Dir')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            for flyStep = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyStep}.ID)
                    break;
                end
                if flyStep == length(allFlyData)
                    allFlyData{numFlies}.ID = flyName;
                    numFlies = numFlies + 1;
                end
            end
        end
    end
end
numFlies = numFlies-1;

% Get the position data, ROI data, and trial ID for all of the flies and all of the trials
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs') & ~strcmpi(fileName(1:3),'PVA') & ~strcmpi(fileName(1:3),'Dir') & ~strcmpi(fileName(end-7:end-4),'Spot')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
            
            % Do some more preprocessing
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 11;
            RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(GROIaveMean,sgolayOrder,sgolayWindow,[],2);
            
            % Match the behavior to the imaging
            num_planes = round(length(positionDat.tFrameGrab)/length(RROIaveMax));
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
            vF = sqrt(diff(OffsetForMatch).^2+diff(OffsetLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));
            vF = sgolayfilt(vF,sgolayOrder,sgolayWindow);
            vSS = diff(OffsetLatMatch)./mean(diff(OffsetRotMatch(:,1)));
            
            positionDatMatch = {};
            positionDatMatch.OffsetRotMatch = OffsetRotMatch;
            positionDatMatch.OffsetForMatch = OffsetForMatch;
            positionDatMatch.OffsetLetMatch = OffsetLatMatch;
            positionDatMatch.vRot = vRot;
            positionDatMatch.vF = vF;
            positionDatMatch.vSS = vSS;
            positionDatMatch.minFG = minFG;
            positionDatMatch.maxFG = maxFG;
            
            for flyID = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyID}.ID)
                    fileParts = strsplit(fileName,'_');
                     if strcmpi(fileParts(5),'1x')
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                    elseif strcmpi(fileParts(5),'2x')
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                    elseif strcmpi(fileParts(5),'Dark')
                        allFlyData{flyID}.dark{str2num(fileName(end-5:end-4))}.positionDatMatch = positionDatMatch;
                        allFlyData{flyID}.dark{str2num(fileName(end-5:end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.dark{str2num(fileName(end-5:end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.dark{str2num(fileName(end-5:end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                        allFlyData{flyID}.dark{str2num(fileName(end-5:end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);  
                    end
                end
            end
        end 
    end
end

%% Make the plots

tStop = 5;
vRThresh = 0.075;
vFThresh = 0.1;
numPts = 17;

WalkStats = figure('units','normalized','outerposition',[0 0 1 1]);

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    vRotAll = [];
    vFAll = [];
    vSSAll =  [];
    pWalkAll = [];
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark)
        vRot = 180/pi*abs(allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot);
        vRotAll = vertcat(vRotAll,vRot);
        vF = 10*allFlyData{flyID}.dark{darkID}.positionDatMatch.vF;
        vFAll = vertcat(vFAll,vF);
        vSSAll = vertcat(vSSAll,10*abs(allFlyData{flyID}.dark{darkID}.positionDatMatch.vSS));
        
        vZero = find(abs(vRot)<vRThresh);
        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);

        vZShort = [];
        for i=1:length(vZBeg)
            if vZEnd(i) - vZBeg(i) < tStop
                vZShort = horzcat(vZShort,[vZBeg(i):vZEnd(i)]);
            end
        end   
        vZero(vZShort) = [];

        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);
        
        vFZero = find(abs(vF)<vFThresh);
        vFZEnd = find(diff(vFZero)>1);
        vFZBeg = vFZEnd+1;
        vFZEnd = vertcat(vFZEnd,length(vFZero));
        vFZBeg = vertcat(1,vFZBeg);

        vFZShort = [];
        for i=1:length(vFZBeg)
            if vFZEnd(i) - vFZBeg(i) < tStop
                vFZShort = horzcat(vFZShort,[vFZBeg(i):vFZEnd(i)]);
            end
        end   
        vFZero(vFZShort) = [];

        vFZEnd = find(diff(vFZero)>1);
        vFZBeg = vFZEnd+1;
        vFZEnd = vertcat(vFZEnd,length(vFZero));
        vFZBeg = vertcat(1,vFZBeg);
        
        
%         figure;
%         hold on;
%         plot(vF,'k');
%         scatter(vFZero,vF(vFZero),'b');
%         scatter(vFZero(vFZEnd),vF(vFZero(vFZEnd)),'r');
%         scatter(vFZero(vFZBeg),vF(vFZero(vFZBeg)),'g');
        
        
        pWalk = 1-union(vZero,vFZero)/length(vRot);
        pWalkAll = vertcat(pWalkAll,pWalk);
    end
    subplot(3,4,1);
    hold on;
    scatter(flyID,mean(vFAll),'k','filled');
    line([flyID flyID],[mean(vFAll)-std(vFAll) mean(vFAll) + std(vFAll)],'Color','k');
    ylim([-1 6]);
    xlabel('fly #');
    ylabel('Forward velocity (mm/sec)');
    xlim([0 numFlies+1]);
    
    subplot(3,4,2);
    hold on;
    scatter(flyID,mean(vSSAll),'k','filled');
    line([flyID flyID],[mean(vSSAll)-std(vSSAll) mean(vSSAll) + std(vSSAll)],'Color','k');
    ylim([-1 4]);
    xlabel('fly #');
    ylabel('Abs sideslip velocity (mm/sec)');
    xlim([0 numFlies+1]);
    
    subplot(3,4,3);
    hold on;
    scatter(flyID,mean(vRotAll),'k','filled');
    line([flyID flyID],[mean(vRotAll)-std(vRotAll) mean(vRotAll) + std(vRotAll)],'Color','k');
    ylim([-20 80]);
    xlabel('fly #');
    ylabel('Abs. rotational velocity (^o/s)');
    xlim([0 numFlies+1]);
    
    subplot(3,4,4);
    hold on;
    scatter(flyID,mean(pWalkAll),'k','filled');
    line([flyID flyID],[mean(pWalkAll)-std(pWalkAll) mean(pWalkAll) + std(pWalkAll)],'Color','k');
    ylim([0 1]);
    xlabel('fly #');
    ylabel('Fraction walking');
    xlim([0 numFlies+1]);
end

set(WalkStats,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(WalkStats,'Results/WalkStats_PB_RT_GW_8135','-dpdf');

%% Save the data to a .csv file
tStop = 5;
vRThresh = 0.075;
vFThresh = 0.1;
numPts = 17;


pWalkStatMean = zeros(numFlies,1);
pWalkStatSD = zeros(numFlies,1);
vRotStatMean = zeros(numFlies,1);
vRotStatSD = zeros(numFlies,1);
vFStatMean = zeros(numFlies,1);
vFStatSD = zeros(numFlies,1);
vSSStatMean = zeros(numFlies,1);
vSSStatSD = zeros(numFlies,1);

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    vRotAll = [];
    vFAll = [];
    vSSAll =  [];
    pWalkAll = [];
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark)
        vRot = 180/pi*abs(allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot);
        vRotAll = vertcat(vRotAll,vRot);
        vF = 10*allFlyData{flyID}.dark{darkID}.positionDatMatch.vF;
        vFAll = vertcat(vFAll,vF);
        vSSAll = vertcat(vSSAll,10*abs(allFlyData{flyID}.dark{darkID}.positionDatMatch.vSS));
        
        vZero = find(abs(vRot)<vRThresh);
        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);

        vZShort = [];
        for i=1:length(vZBeg)
            if vZEnd(i) - vZBeg(i) < tStop
                vZShort = horzcat(vZShort,[vZBeg(i):vZEnd(i)]);
            end
        end   
        vZero(vZShort) = [];

        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);
        
        vFZero = find(abs(vF)<vFThresh);
        vFZEnd = find(diff(vFZero)>1);
        vFZBeg = vFZEnd+1;
        vFZEnd = vertcat(vFZEnd,length(vFZero));
        vFZBeg = vertcat(1,vFZBeg);

        vFZShort = [];
        for i=1:length(vFZBeg)
            if vFZEnd(i) - vFZBeg(i) < tStop
                vFZShort = horzcat(vFZShort,[vFZBeg(i):vFZEnd(i)]);
            end
        end   
        vFZero(vFZShort) = [];

        vFZEnd = find(diff(vFZero)>1);
        vFZBeg = vFZEnd+1;
        vFZEnd = vertcat(vFZEnd,length(vFZero));
        vFZBeg = vertcat(1,vFZBeg);
        
        
%         figure;
%         hold on;
%         plot(vF,'k');
%         scatter(vFZero,vF(vFZero),'b');
%         scatter(vFZero(vFZEnd),vF(vFZero(vFZEnd)),'r');
%         scatter(vFZero(vFZBeg),vF(vFZero(vFZBeg)),'g');
        
        
        pWalk = 1-union(vZero,vFZero)/length(vRot);
        pWalkAll = vertcat(pWalkAll,pWalk);
    end
    
    pWalkStatMean(flyID) = mean(pWalkAll);
    pWalkStatSD(flyID) = std(pWalkAll);
    vRotStatMean(flyID) = mean(vRotAll);
    vRotStatSD(flyID) = std(vRotAll);
    vFStatMean(flyID) = mean(vFAll);
    vFStatSD(flyID) =  std(vFAll);
    vSSStatMean(flyID) = mean(vSSAll);
    vSSStatSD(flyID) = std(vSSAll);
end


pms = {};
for i = 1:length(pWalkStatMean)
    pms{i} = '\pm';
end
pms = char(pms);

commas = {};
for i = 1:length(pWalkStatMean)
    commas{i} = ',';
end
commas = char(commas);

AllDatCSV = horzcat(...
    num2str([1:numFlies]'),commas,...
    horzcat(num2str(pWalkStatMean,2),pms,num2str(pWalkStatSD,2)),commas,...
    horzcat(num2str(vRotStatMean,2),pms,num2str(vRotStatSD,2)),commas,...
    horzcat(num2str(vFStatMean,2),pms,num2str(vFStatSD,2)),commas,...
    horzcat(num2str(vSSStatMean,2),pms,num2str(vSSStatSD,2)));

setName = strcat(region,'_G',greenLine,'_R',redLine);
dlmwrite(strcat('D:\Imaging\2Color\Stats\',setName,'.csv'),AllDatCSV,'delimiter','');

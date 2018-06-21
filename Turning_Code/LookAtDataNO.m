%% Clear out old data
clear;
clc;

%% Load all of the data
% Specify directories
dirRoot = 'D:\Imaging\2Color';

imDir{1} = '20160423-Flip';
imDir{2} = '20160425-Flip';
imDir{3} = '20160426-Flip';
imDir{4} = '20160429';
imDir{5} = '20160505';
redLine = '60D05';
greenLine = '37F06';

cd(strcat(dirRoot,'\',imDir{1}));

% Get the fly IDs
allFlyData{1}.ID = 'Empty';
numFlies = 1;
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-5:end),'NO.mat')
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
        if strcmpi(fileName(end-5:end),'NO.mat')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
            
            % Do some more preprocessing
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 11;
            RROIaveFilt = sgolayfilt(RROIave,sgolayOrder,sgolayWindow,[],2);
            GROIaveFilt = sgolayfilt(GROIave,sgolayOrder,sgolayWindow,[],2);
            
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
            vF = sqrt(diff(OffsetForMatch).^2+diff(OffsetLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));
            vF = sgolayfilt(vF,sgolayOrder,sgolayWindow);
            
            positionDatMatch = {};
            positionDatMatch.OffsetRotMatch = OffsetRotMatch;
            positionDatMatch.OffsetForMatch = OffsetForMatch;
            positionDatMatch.OffsetLetMatch = OffsetLatMatch;
            positionDatMatch.vRot = vRot;
            positionDatMatch.vF = vF;
            positionDatMatch.minFG = minFG;
            positionDatMatch.maxFG = maxFG;
            
            for flyID = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyID}.ID)
                    fileParts = strsplit(fileName,'_');
                         if strcmpi(fileParts(5),'1x')
                            allFlyData{flyID}.gainOne{str2num(fileName(end-7))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.gainOne{str2num(fileName(end-7))}.RROIave = RROIaveFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainOne{str2num(fileName(end-7))}.GROIave = GROIaveFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'2x')
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-7))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-7))}.RROIave = RROIaveFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-7))}.GROIave = GROIaveFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'Dark')
                            allFlyData{flyID}.dark{str2num(fileName(end-7))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.dark{str2num(fileName(end-7))}.RROIave = RROIaveFilt(:,minFG:maxFG);
                            allFlyData{flyID}.dark{str2num(fileName(end-7))}.GROIave = GROIaveFilt(:,minFG:maxFG);
                    end
                end
            end
        end 
    end
end

%% Create a scatter plot of the correlations

CCs = {};

offset = 3;

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        ROISig = allFlyData{flyID}.dark{darkID}.GROIave(1,:)-allFlyData{flyID}.dark{darkID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2+offset:end),allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot(1:end-offset));
        CCs{flyID}.dark(darkID) = CCNow(2,1);
    end
    
    for gainOneID = 1:length(allFlyData{flyID}.gainOne);
        ROISig = allFlyData{flyID}.gainOne{gainOneID}.GROIave(1,:)-allFlyData{flyID}.gainOne{gainOneID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2+offset:end),allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.vRot(1:end-offset));
        CCs{flyID}.gainOne(gainOneID) = CCNow(2,1);
    end
    
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo);
        ROISig = allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(1,:)-allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2+offset:end),allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.vRot(1:end-offset));
        CCs{flyID}.gainTwo(gainTwoID) = CCNow(2,1);
    end
end

NOCCs = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.dark),1)*[1 0 0];
    scatter(zeros(length(CCs{flyID}.dark),1)+flyID,CCs{flyID}.dark,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.dark) median(CCs{flyID}.dark)], 'Color','r','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
ylabel('correlation (\DeltaF/F vs. v_{R})');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('dark');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,2);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainOne),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainOne),1)+flyID,CCs{flyID}.gainOne,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainOne) median(CCs{flyID}.gainOne)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (1x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,3);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainTwo),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainTwo),1)+flyID,CCs{flyID}.gainTwo,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainTwo) median(CCs{flyID}.gainTwo)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (2x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


set(NOCCs,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(NOCCs,'..\Results\NOCCs','-dpdf');

%% Create a scatter plot of the correlations

CCs = {};

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        ROISig = allFlyData{flyID}.dark{darkID}.GROIave(1,:)-allFlyData{flyID}.dark{darkID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2:end),allFlyData{flyID}.dark{darkID}.positionDatMatch.vF);
        CCs{flyID}.dark(darkID) = CCNow(2,1);
    end
    
    for gainOneID = 1:length(allFlyData{flyID}.gainOne);
        ROISig = allFlyData{flyID}.gainOne{gainOneID}.GROIave(1,:)-allFlyData{flyID}.gainOne{gainOneID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2:end),allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.vF);
        CCs{flyID}.gainOne(gainOneID) = CCNow(2,1);
    end
    
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo);
        ROISig = allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(1,:)-allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(2,:);
        CCNow = corrcoef(ROISig(2:end),allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.vF);
        CCs{flyID}.gainTwo(gainTwoID) = CCNow(2,1);
    end
end

NOCCs = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.dark),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.dark),1)+flyID,CCs{flyID}.dark,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.dark) median(CCs{flyID}.dark)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
ylabel('correlation (\DeltaF/F vs. v_{F})');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('dark');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,2);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainOne),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainOne),1)+flyID,CCs{flyID}.gainOne,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainOne) median(CCs{flyID}.gainOne)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (1x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,3);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainTwo),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainTwo),1)+flyID,CCs{flyID}.gainTwo,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainTwo) median(CCs{flyID}.gainTwo)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (2x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


set(NOCCs,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(NOCCs,'..\Results\NOCCs_vF','-dpdf');

%% Mean, SD, median, quartile ranges for CCs
allDarkCCs =[];
for fly=1:length(CCs)
    allDarkCCs = horzcat(allDarkCCs,CCs{fly}.dark);
end

mean(allDarkCCs)
std(allDarkCCs)


quantile(allDarkCCs,0.25)
median(allDarkCCs)
quantile(allDarkCCs,0.75)

%% Create a scatter plot of the correlations - convoluting the vR signal
% GtOn = 0.42/log(2);
% GtOff = 0.33/log(2);

GtOn = 0.092/log(2);
GtOff = 0.44/log(2);
    
CCs = {};

offset = 0;

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        ROISig = allFlyData{flyID}.dark{darkID}.GROIave(1,:)-allFlyData{flyID}.dark{darkID}.GROIave(2,:);
        
        vRot = allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot(1:end);
        tAll = allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
%         
%         vRotCCW = vRot;
%         vRotCCW(find(vRotCCW<0)) = 0;
%         vRotCW = -vRot;
%         vRotCW(find(vRotCW<0)) = 0;
%         
%         vRotCCWGoingUp = find(diff(vRotCCW)>0);
%         vRotCCWGoingDown = find(diff(vRotCCW)<0);
%         vRotCCWUp = zeros(length(vRotCCW),1);
%         vRotCCWUp(vRotCCWGoingUp)= vRotCCW(vRotCCWGoingUp);
%         vRotCCWDown = zeros(length(vRotCCW),1);
%         vRotCCWDown(vRotCCWGoingDown)= vRotCCW(vRotCCWGoingDown);
%         vRotCCWConv = conv(vRotCCWUp,exp(-tAll/GtOn)) + conv(vRotCCWDown,exp(-tAll/GtOff));
%         
%         vRotCWGoingUp = find(diff(vRotCW)>0);
%         vRotCWGoingDown = find(diff(vRotCW)<0);
%         vRotCWUp = zeros(length(vRotCW),1);
%         vRotCWUp(vRotCWGoingUp)= vRotCW(vRotCWGoingUp);
%         vRotCWDown = zeros(length(vRotCW),1);
%         vRotCWDown(vRotCWGoingDown)= vRotCW(vRotCWGoingDown);
%         vRotCWConv = conv(vRotCWUp,exp(-tAll/GtOn)) + conv(vRotCWDown,exp(-tAll/GtOff));
%         
%         vRotConv = vRotCCWConv-vRotCWConv;
        vRotConv = conv(vRot,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        
        
        [CCNow, R] = corrcoef(ROISig((2+offset):end),vRotConv(1:(length(vRot)-offset)),'alpha',0.0000001);
        CCs{flyID}.dark(darkID) = CCNow(2,1);
    end
    
    % Plot parameters for the dark data
    for gainOneID = 1:length(allFlyData{flyID}.gainOne);
        ROISig = allFlyData{flyID}.gainOne{gainOneID}.GROIave(1,:)-allFlyData{flyID}.gainOne{gainOneID}.GROIave(2,:);
        
        vRot = allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.vRot(1:end);
        tAll = allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
%         
%         vRotCCW = vRot;
%         vRotCCW(find(vRotCCW<0)) = 0;
%         vRotCW = -vRot;
%         vRotCW(vRotCW<0) = 0;
%         
%         vRotCCWGoingUp = find(diff(vRotCCW)>0);
%         vRotCCWGoingDown = find(diff(vRotCCW)<0);
%         vRotCCWUp = zeros(length(vRotCCW),1);
%         vRotCCWUp(vRotCCWGoingUp)= vRotCCW(vRotCCWGoingUp);
%         vRotCCWDown = zeros(length(vRotCCW),1);
%         vRotCCWDown(vRotCCWGoingDown)= vRotCCW(vRotCCWGoingDown);
%         vRotCCWConv = conv(vRotCCWUp,exp(-tAll/GtOn)) + conv(vRotCCWDown,exp(-tAll/GtOff));
%         
%         vRotCWGoingUp = find(diff(vRotCW)>0);
%         vRotCWGoingDown = find(diff(vRotCW)<0);
%         vRotCWUp = zeros(length(vRotCW),1);
%         vRotCWUp(vRotCWGoingUp)= vRotCW(vRotCWGoingUp);
%         vRotCWDown = zeros(length(vRotCW),1);
%         vRotCWDown(vRotCWGoingDown)= vRotCW(vRotCWGoingDown);
%         vRotCWConv = conv(vRotCWUp,exp(-tAll/GtOn)) + conv(vRotCWDown,exp(-tAll/GtOff));
%         
%         vRotConv = vRotCCWConv-vRotCWConv;
        vRotConv = conv(vRot,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        
        
        CCNow = corrcoef(ROISig((2+offset):end),vRotConv(1:(length(vRot)-offset)));
        CCs{flyID}.gainOne(gainOneID) = CCNow(2,1);
    end
    
    % Plot parameters for the dark data
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo);
        ROISig = allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(1,:)-allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(2,:);
        
        vRot = allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.vRot(1:end);
        tAll = allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
%         
%         vRotCCW = vRot;
%         vRotCCW(find(vRotCCW<0)) = 0;
%         vRotCW = -vRot;
%         vRotCW(find(vRotCW<0)) = 0;
%         
%         vRotCCWGoingUp = find(diff(vRotCCW)>0);
%         vRotCCWGoingDown = find(diff(vRotCCW)<0);
%         vRotCCWUp = zeros(length(vRotCCW),1);
%         vRotCCWUp(vRotCCWGoingUp)= vRotCCW(vRotCCWGoingUp);
%         vRotCCWDown = zeros(length(vRotCCW),1);
%         vRotCCWDown(vRotCCWGoingDown)= vRotCCW(vRotCCWGoingDown);
%         vRotCCWConv = conv(vRotCCWUp,exp(-tAll/GtOn)) + conv(vRotCCWDown,exp(-tAll/GtOff));
%         
%         vRotCWGoingUp = find(diff(vRotCW)>0);
%         vRotCWGoingDown = find(diff(vRotCW)<0);
%         vRotCWUp = zeros(length(vRotCW),1);
%         vRotCWUp(vRotCWGoingUp)= vRotCW(vRotCWGoingUp);
%         vRotCWDown = zeros(length(vRotCW),1);
%         vRotCWDown(vRotCWGoingDown)= vRotCW(vRotCWGoingDown);
%         vRotCWConv = conv(vRotCWUp,exp(-tAll/GtOn)) + conv(vRotCWDown,exp(-tAll/GtOff));
%         
%         vRotConv = vRotCCWConv-vRotCWConv;
        vRotConv = conv(vRot,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        
        
        CCNow = corrcoef(ROISig((2+offset):end),vRotConv(1:(length(vRot)-offset)));
        CCs{flyID}.gainTwo(gainTwoID) = CCNow(2,1);
    end
end

NOCCs = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.dark),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.dark),1)+flyID,CCs{flyID}.dark,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [mean(CCs{flyID}.dark) mean(CCs{flyID}.dark)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
ylabel('correlation (\DeltaF/F vs. v_{R})');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('dark');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,2);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainOne),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainOne),1)+flyID,CCs{flyID}.gainOne,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [mean(CCs{flyID}.gainOne) mean(CCs{flyID}.gainOne)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (1x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,3);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainTwo),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainTwo),1)+flyID,CCs{flyID}.gainTwo,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [mean(CCs{flyID}.gainTwo) mean(CCs{flyID}.gainTwo)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (2x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');

% set(NOCCs,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(NOCCs,'..\Results\NOCCs_Convolved','-dpdf');

%% Create a scatter plot of the correlations - convoluting the vF signal
% GtOn = 0.42/log(2);
% GtOff = 0.33/log(2);

GtOn = 0.092/log(2);
GtOff = 0.44/log(2);

CCs = {};

offset = 0;

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        ROISig = allFlyData{flyID}.dark{darkID}.GROIave(1,:)+allFlyData{flyID}.dark{darkID}.GROIave(2,:);
        
        vF = allFlyData{flyID}.dark{darkID}.positionDatMatch.vF;
        tAll = allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
%         
%         vFGoingUp = find(diff(vF)>0);
%         vFGoingDown = find(diff(vF)<0);
%         vFUp = zeros(length(vF),1);
%         vFUp(vFGoingUp)= vF(vFGoingUp);
%         vFDown = zeros(length(vF),1);
%         vFDown(vFGoingDown)= vF(vFGoingDown);
%         vFConv = conv(vFUp,exp(-tAll/GtOn)) + conv(vFDown,exp(-tAll/GtOff));
        vFConv = conv(vF,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        
        CCNow = corrcoef(ROISig((2+offset):end),vFConv(1:(length(vF)-offset)));
        CCs{flyID}.dark(darkID) = CCNow(2,1);
    end
    
    % Plot parameters for the dark data
    for gainOneID = 1:length(allFlyData{flyID}.gainOne);
        ROISig = allFlyData{flyID}.gainOne{gainOneID}.GROIave(1,:)-allFlyData{flyID}.gainOne{gainOneID}.GROIave(2,:);
        
        vF = allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.vF;
        tAll = allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
%         
%        vFGoingUp = find(diff(vF)>0);
%         vFGoingDown = find(diff(vF)<0);
%         vFUp = zeros(length(vF),1);
%         vFUp(vFGoingUp)= vF(vFGoingUp);
%         vFDown = zeros(length(vF),1);
%         vFDown(vFGoingDown)= vF(vFGoingDown);
%         vFConv = conv(vFUp,exp(-tAll/GtOn)) + conv(vFDown,exp(-tAll/GtOff));
        vFConv = conv(vF,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        
        CCNow = corrcoef(ROISig((2+offset):end),vFConv(1:(length(vF)-offset)));    
        CCs{flyID}.gainOne(gainOneID) = CCNow(2,1);
    end
    
    % Plot parameters for the dark data
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo);
        ROISig = allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(1,:)-allFlyData{flyID}.gainTwo{gainTwoID}.GROIave(2,:);
        
        vF = allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.vF;
        tAll = allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(1:end-offset,1);
        
%         vFGoingUp = find(diff(vF)>0);
%         vFGoingDown = find(diff(vF)<0);
%         vFUp = zeros(length(vF),1);
%         vFUp(vFGoingUp)= vF(vFGoingUp);
%         vFDown = zeros(length(vF),1);
%         vFDown(vFGoingDown)= vF(vFGoingDown);
%         vFConv = conv(vFUp,exp(-tAll/GtOn)) + conv(vFDown,exp(-tAll/GtOff));
        vFConv = conv(vF,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));        
        CCNow = corrcoef(ROISig((2+offset):end),vFConv(1:(length(vF)-offset)));
        CCs{flyID}.gainTwo(gainTwoID) = CCNow(2,1);
    end
end

NOCCs = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.dark),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.dark),1)+flyID,CCs{flyID}.dark,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [mean(CCs{flyID}.dark) mean(CCs{flyID}.dark)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
ylabel('correlation (\DeltaF/F vs. v_{R})');
xlim([ 0 (numFlies+1)]);
ylim([0 1]);
title('dark');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,2);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainOne),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainOne),1)+flyID,CCs{flyID}.gainOne,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainOne) median(CCs{flyID}.gainOne)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (1x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');


subplot(2,3,3);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.gainTwo),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.gainTwo),1)+flyID,CCs{flyID}.gainTwo,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.gainTwo) median(CCs{flyID}.gainTwo)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('stripe (2x gain)');
line([0 (numFlies+1)], [0 0],'Color','k');

% set(NOCCs,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(NOCCs,'..\Results\NOCCs_VF_Convolved','-dpdf');

%% Plot an example trace
flyID = 8;
trial = 5;

% GtOn = 0.42/log(2);
% GtOff = 0.33/log(2);
GtOn = 0.092/log(2);
GtOff = 0.44/log(2);
    

figure;
subplot(3,1,1);
hold on;
plot(allFlyData{flyID}.dark{trial}.GROIave(1,2:end));
plot(allFlyData{flyID}.dark{trial}.GROIave(2,2:end));
line([0 length(allFlyData{flyID}.dark{trial}.GROIave(1,2:end))],[ 1 1]);

vRot = allFlyData{flyID}.dark{trial}.positionDatMatch.vRot(1:end);
tAll = allFlyData{flyID}.dark{trial}.positionDatMatch.OffsetRotMatch(:,1);

% vRotCCW = vRot;
% vRotCCW(find(vRotCCW<0)) = 0;
% vRotCW = -vRot;
% vRotCW(find(vRotCW<0)) = 0;
% 
% vRotCCWGoingUp = find(diff(vRotCCW)>0);
% vRotCCWGoingDown = find(diff(vRotCCW)<0);
% vRotCCWUp = zeros(length(vRotCCW),1);
% vRotCCWUp(vRotCCWGoingUp)= vRotCCW(vRotCCWGoingUp);
% vRotCCWDown = zeros(length(vRotCCW),1);
% vRotCCWDown(vRotCCWGoingDown)= vRotCCW(vRotCCWGoingDown);
% vRotCCWConv = conv(vRotCCWUp,exp(-tAll/GtOn)) + conv(vRotCCWDown,exp(-tAll/GtOff));
% 
% vRotCWGoingUp = find(diff(vRotCW)>0);
% vRotCWGoingDown = find(diff(vRotCW)<0);
% vRotCWUp = zeros(length(vRotCW),1);
% vRotCWUp(vRotCWGoingUp)= vRotCW(vRotCWGoingUp);
% vRotCWDown = zeros(length(vRotCW),1);
% vRotCWDown(vRotCWGoingDown)= vRotCW(vRotCWGoingDown);
% vRotCWConv = conv(vRotCWUp,exp(-tAll/GtOn)) + conv(vRotCWDown,exp(-tAll/GtOff));
% 
% vRotConv = vRotCCWConv-vRotCWConv;
vRotConv = conv(vRot,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));

subplot(3,1,2);
plot(vRotConv(1:length(vRotConv)/2));

subplot(3,1,3);
plot(-allFlyData{flyID}.dark{trial}.GROIave(2,2:end)+allFlyData{flyID}.dark{trial}.GROIave(1,2:end));
line([0 length(allFlyData{flyID}.dark{trial}.GROIave(1,2:end))],[0 0]);
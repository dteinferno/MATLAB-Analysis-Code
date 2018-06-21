%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to analyze all of the 2Color PB data (already processed with ProcessAllData2ColorPB) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear out old data
clear;
clc;

%% Load all of the data
% Specify directories
dirRoot = 'D:\Imaging\2Color';
% dirRoot = 'D:\Imaging\2Color\20160718';
% dirRoot = 'D:\Imaging\2Color\NegativeGain\';

% imDir{1} = '20170118';
% redLine = '60D05';
% greenLine = '37F06';


% imDir{1} = '20160617';
% redLine = '60D05';
% greenLine = 'VT20739';
% 
% imDir{1} = '20160621';
% imDir{2} = '20160908';
% redLine = '37F06';
% redSpan = [1:8 11:18];
% greenLine = '60D05';
% greenSpan = [2:17];
% 
imDir{1} = '20160622';
imDir{2} = '20160907';
imDir{3} = '20160908-Flip';
redLine = '60D05';
redSpan = [2:17];
greenLine = '37F06';
greenSpan = [1:8 11:18];

% imDir{1} = '20160627';
% imDir{2} = '20160630';
% redLine = 'VT20739';
% greenLine = '60D05';

% imDir{1} = 'VTRed';
% redLine = 'VT20739';
% greenLine = '60D05';

% imDir{1} = 'PEN2\20160718';
% redLine = '60D05';
% redSpan = [2:17];
% greenLine = 'VT20739';
% greenSpan = [1:8 11:18];

% imDir{1} = '20160719';
% imDir{2} = '20160723';
% imDir{3} = '20160725';
% redLine = 'VT45661';
% greenLine = '60D05';

% imDir{1} = '20160903';
% redLine = 'VT8135';
% redSpan = [1:8 11:18];
% greenLine = '60D05';
% greenSpan = [2:17];

% dirRoot = 'D:\Imaging\37F06\20160928';
% imDir{1} = 'DMD';
% redLine = '60D05';
% redSpan = [2:17];
% greenLine = '37F06';
% greenSpan = [1:8 11:18];

cd(strcat(dirRoot,'\',imDir{1}));

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
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                    elseif strcmpi(fileParts(5),'-1x')
                        allFlyData{flyID}.gainNOne{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                        allFlyData{flyID}.gainNOne{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainNOne{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainNOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
                        allFlyData{flyID}.gainNOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
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

%% Plot a summary for individual dark, 1x, and 2x trials

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID});

            % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(allFlyData{flyID}.dark{darkID}.RROIaveMax);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(allFlyData{flyID}.dark{darkID}.GROIaveMax);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(allFlyData{flyID}.dark{darkID}.RROIaveMax) 3]);
            overlayIm(:,:,1) = (allFlyData{flyID}.dark{darkID}.RROIaveMax-...
                min(min(allFlyData{flyID}.dark{darkID}.RROIaveMax)))./...
                (max(max(allFlyData{flyID}.dark{darkID}.RROIaveMax))-min(min(allFlyData{flyID}.dark{darkID}.RROIaveMax)));
            overlayIm(:,:,2) = (allFlyData{flyID}.dark{darkID}.GROIaveMax-...
                min(min(allFlyData{flyID}.dark{darkID}.GROIaveMax)))./...
                (max(max(allFlyData{flyID}.dark{darkID}.GROIaveMax))-min(min(allFlyData{flyID}.dark{darkID}.GROIaveMax)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            subplot(4,1,4);
            plot(allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(:,1),...
                allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(:,2),'k','LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(1,1) allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(end,1)]);
            ylim([-pi pi]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_Dark_Summary_',num2str(darkID)),'-dpdf');
            print(RGcomp,strcat(allFlyData{flyID}.ID,'_Dark_Summary_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
    
    % Plot parameters for the 1x gain data
    for gainOneID = 1:length(allFlyData{flyID}.gainOne)
        if ~isempty(allFlyData{flyID}.gainOne{gainOneID});
             % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax) 3]);
            overlayIm(:,:,1) = (allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax-...
                min(min(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax)))./...
                (max(max(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax))-min(min(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax)));
            overlayIm(:,:,2) = (allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax-...
                min(min(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax)))./...
                (max(max(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax))-min(min(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            subplot(4,1,4);
            plot(allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(:,1),...
                allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(:,2),'LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(1,1) allFlyData{flyID}.gainOne{gainOneID}.positionDatMatch.OffsetRotMatch(end,1)]);
            ylim([-pi pi]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_gainOne_Summary_',num2str(gainOneID)),'-dpdf');
            print(RGcomp,strcat(allFlyData{flyID}.ID,'_gainOne_Summary_',num2str(gainOneID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
    
    % Plot parameters for the 2x gain data
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo)
        if ~isempty(allFlyData{flyID}.gainTwo{gainTwoID});
             % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax) 3]);
            overlayIm(:,:,1) = (allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax-...
                min(min(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax)))./...
                (max(max(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax))-min(min(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax)));
            overlayIm(:,:,2) = (allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax-...
                min(min(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax)))./...
                (max(max(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax))-min(min(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            stripeMatch = mod(2*UnWrap(allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(:,2),2,0),2*pi)-pi;
            
            subplot(4,1,4);
            hold on;
            plot(allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(:,1),...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(:,2),'LineWidth',2,'LineStyle','--');
            plot(allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(:,1),...
                stripeMatch,'LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(1,1) allFlyData{flyID}.gainTwo{gainTwoID}.positionDatMatch.OffsetRotMatch(end,1)]);
            ylim([-pi pi]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_gainTwo_Summary_',num2str(gainTwoID)),'-dpdf');
            print(RGcomp,strcat(allFlyData{flyID}.ID,'_gainTwo_Summary_',num2str(gainTwoID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

%% Plot the peak maximum vs. the rotational velocity

% Specify the range of rotational velocities to consider
vRMin = 0;
vRMax = 720;
vRSpan = 30;
[LPB.dark.dataR, LPB.dark.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(1:8),greenSpan(1:8), 'dark');
[LPB.gainOne.dataR, LPB.gainOne.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(1:8),greenSpan(1:8), 'gainOne');
[LPB.gainTwo.dataR, LPB.gainTwo.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(1:8),greenSpan(1:8), 'gainTwo');
[RPB.dark.dataR, RPB.dark.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(9:16),greenSpan(9:16), 'dark');
[RPB.gainOne.dataR, RPB.gainOne.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(9:16),greenSpan(9:16), 'gainOne');
[RPB.gainTwo.dataR, RPB.gainTwo.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(9:16),greenSpan(9:16), 'gainTwo');

plotMax = 180;
% Look at the red channel   
vRotAmpR = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    % Look at the left PB in the dark
    subplot(6,numFlies,flyID);
    hold on;
    TurnPlot(LPB.dark.dataR{flyID},vRMin,vRMax,vRSpan, plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    title(strcat(redLine,'-Fly ',num2str(flyID)));
    text(vRMax*pi/180*0.8,1.8,'CW','Color','c','FontSize',14);
    text(vRMax*pi/180*0.8,1.6,'CCW','Color','m','FontSize',14);
    ylabel('Left PB, dark');
    
    % Look at the right PB in the dark
    subplot(6,numFlies,flyID+numFlies);
    hold on;
    TurnPlot(RPB.dark.dataR{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    ylabel('Right PB, dark');
    
    % Look at the left PB with a gain=1 stripe
    subplot(6,numFlies,flyID+numFlies*2);
    hold on;
    TurnPlot(LPB.gainOne.dataR{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    ylabel('Left PB, 1x gain');
    
    % Look at the right PB with a gain=1 stripe
    subplot(6,numFlies,flyID+numFlies*3);
    hold on;
    TurnPlot(RPB.gainOne.dataR{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    ylabel('Right PB, 1x gain');
    
    % Look at the left PB with a gain=2 stripe
    subplot(6,numFlies,flyID+numFlies*4);
    hold on;
    TurnPlot(LPB.gainTwo.dataR{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    ylabel('Left PB, 2x gain');
    
    % Look at the right PB with a gain=2 stripe
    subplot(6,numFlies,flyID+numFlies*5);
    hold on;
    TurnPlot(RPB.gainTwo.dataR{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 2]);
    xlabel('Speed (rad/sec)');
    ylabel('Right PB, 2x gain');
end

% Look at the green channel   
vRotAmpG = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    % Look at the left PB in the dark
    subplot(6,numFlies,flyID);
    hold on;
    TurnPlot(LPB.dark.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    title(strcat(greenLine,'-Fly ',num2str(flyID)));
    text(vRMax*pi/180*0.8,2.8,'CW','Color','c','FontSize',14);
    text(vRMax*pi/180*0.8,2.4,'CCW','Color','m','FontSize',14);
    ylabel('Left PB, dark');
    
    % Look at the right PB in the dark
    subplot(6,numFlies,flyID+numFlies);
    hold on;
    TurnPlot(RPB.dark.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    ylabel('Right PB, dark');
    
    % Look at the left PB with a gain=1 stripe
    subplot(6,numFlies,flyID+numFlies*2);
    hold on;
    TurnPlot(LPB.gainOne.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    ylabel('Left PB, 1x gain');
    
    % Look at the right PB with a gain=1 stripe
    subplot(6,numFlies,flyID+numFlies*3);
    hold on;
    TurnPlot(RPB.gainOne.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    ylabel('Right PB, 1x gain');
    
    % Look at the left PB with a gain=2 stripe
    subplot(6,numFlies,flyID+numFlies*4);
    hold on;
    TurnPlot(LPB.gainTwo.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    ylabel('Left PB, 2x gain');
    
    % Look at the right PB with a gain=2 stripe
    subplot(6,numFlies,flyID+numFlies*5);
    hold on;
    TurnPlot(RPB.gainTwo.dataG{flyID},vRMin,vRMax,vRSpan,plotMax/vRSpan);
    xlim([0 plotMax*pi/180]);
    ylim([1 3]);
    xlabel('Speed (rad/sec)');
    ylabel('Right PB, 2x gain');
end

% set(vRotAmpR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(vRotAmpR,strcat(redLine,'_vRotVsMax'),'-dpdf');
% set(vRotAmpG,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(vRotAmpG,strcat(greenLine,'_vRotVsMax'),'-dpdf');

%% Align the peaks and plot the normalized mean values as sorted by velocity
glomShift = 5;
vRMin = 0;
vRMax = 720;
vRSpan = 30;
[LPB.dark.dataR, RPB.dark.dataR, LPB.dark.dataG, RPB.dark.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'dark',glomShift);
[LPB.gainOne.dataR, RPB.gainOne.dataR, LPB.gainOne.dataG, RPB.gainOne.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainOne',glomShift);
[LPB.gainTwo.dataR, RPB.gainTwo.dataR, LPB.gainTwo.dataG, RPB.gainTwo.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainTwo',glomShift);

vRBinNum = round((vRMax-vRMin)/vRSpan);

PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

plotRange = 180/vRSpan;
for flyID = 1:length(allFlyData)
    
    LPBDataRAllStop = horzcat(horzcat(LPB.dark.dataR{flyID}.Stop,LPB.gainOne.dataR{flyID}.Stop),LPB.gainTwo.dataR{flyID}.Stop);
    RPBDataRAllStop = horzcat(horzcat(RPB.dark.dataR{flyID}.Stop,RPB.gainOne.dataR{flyID}.Stop),RPB.gainTwo.dataR{flyID}.Stop);
    LPBDataGAllStop = horzcat(horzcat(LPB.dark.dataG{flyID}.Stop,LPB.gainOne.dataG{flyID}.Stop),LPB.gainTwo.dataG{flyID}.Stop);
    RPBDataGAllStop = horzcat(horzcat(RPB.dark.dataG{flyID}.Stop,RPB.gainOne.dataG{flyID}.Stop),RPB.gainTwo.dataG{flyID}.Stop);
    
    % Plot the dark data
    LRStop = mean(LPBDataRAllStop,2);
    LRStop = (LRStop-min(LRStop(redSpan(1:8))))./(max(LRStop)-min(LRStop(redSpan(1:8))));
    RRStop = mean(RPBDataRAllStop,2);
    RRStop = (RRStop-min(RRStop(greenSpan(1:8))))./(max(RRStop)-min(RRStop(greenSpan(1:8))));
    RStop = vertcat(LRStop,RRStop);
    
    LGStop = mean(LPBDataGAllStop,2);
    LGStop = (LGStop-min(LGStop(greenSpan(1:8))))./(max(LGStop)-min(LGStop(greenSpan(1:8))));
    RGStop = median(RPBDataGAllStop,2);
    RGStop = (RGStop-min(RGStop(redSpan(1:8))))./(max(RGStop)-min(RGStop(redSpan(1:8))));
    GStop = vertcat(LGStop,RGStop);
    
    LRPk = find(LRStop == max(LRStop));
    LGPk = find(LGStop == max(LGStop));
    pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
    RRPk = find(RRStop == max(RRStop));
    RGPk = find(RGStop == max(RGStop));
    pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
    
    subplot(4*numFlies,plotRange+1,1+4*(plotRange+1)*(flyID-1))
    hold on;
    LBin = mean(LPBDataRAllStop,2);
    RBin = mean(RPBDataRAllStop,2);
    imagesc(vertcat(LBin,RBin)'-1);   
    line([(redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
        (redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    line([(11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
        (11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    text(2,0,num2str(pkDiffL));
    text(11,0,num2str(pkDiffR));
    title(strcat('Left PB, Fly ',num2str(flyID)));
    axis off;
    caxis([0 0.5]);
    
    subplot(4*numFlies,plotRange+1,1+4*(plotRange+1)*(flyID-1)+plotRange+1)
    LBin =  mean(LPBDataGAllStop,2);
    RBin = mean(RPBDataGAllStop,2);
    imagesc(vertcat(LBin,RBin)'-1);
    line([(greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
        (greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    line([(11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
        (11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
        [0.5 1.5],'Color','k','LineWidth',2);
    axis off;
    caxis([0 0.5]);
    
    
    for binID=1:plotRange
        LPBDataRAllCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CW{binID},LPB.gainOne.dataR{flyID}.CW{binID}),LPB.gainTwo.dataR{flyID}.CW{binID});
        RPBDataRAllCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CW{binID},RPB.gainOne.dataR{flyID}.CW{binID}),RPB.gainTwo.dataR{flyID}.CW{binID});
        LPBDataGAllCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CW{binID},LPB.gainOne.dataG{flyID}.CW{binID}),LPB.gainTwo.dataG{flyID}.CW{binID});
        RPBDataGAllCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CW{binID},RPB.gainOne.dataG{flyID}.CW{binID}),RPB.gainTwo.dataG{flyID}.CW{binID});
        
        if ~isempty(LPBDataRAllCW)
            LRCW = mean(LPBDataRAllCW,2);
            LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
            RRCW = mean(RPBDataRAllCW,2);
            RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
            RCW = vertcat(LRCW,RRCW);
            LGCW = mean(LPBDataGAllCW,2);
            LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
            RGCW = mean(RPBDataGAllCW,2);
            RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
            GCW = vertcat(LGCW,RGCW);

            LRPk = find(LRCW == max(LRCW));
            LGPk = find(LGCW == max(LGCW));
            pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
            RRPk = find(RRCW == max(RRCW));
            RGPk = find(RGCW == max(RGCW));
            pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
        
        
            subplot(4*numFlies,plotRange+1,1+binID+4*(plotRange+1)*(flyID-1))
            LBin = mean(LPBDataRAllCW,2);
            RBin = mean(RPBDataRAllCW,2);
            imagesc(vertcat(LBin,RBin)'-1);
            line([(redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
                (redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            line([(11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
                (11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            text(2,0,num2str(pkDiffL));
            text(11,0,num2str(pkDiffR));
            axis off;
            caxis([0 0.5]);

            subplot(4*numFlies,plotRange+1,1+binID+4*(plotRange+1)*(flyID-1)+plotRange+1)
            LBin = mean(LPBDataGAllCW,2);
            RBin = mean(RPBDataGAllCW,2);
            imagesc(vertcat(LBin,RBin)'-1);
            line([(greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
                (greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            line([(11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
                (11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            axis off;
            caxis([0 0.5]);
        end
        
        LPBDataRAllCCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CCW{binID},LPB.gainOne.dataR{flyID}.CCW{binID}),LPB.gainTwo.dataR{flyID}.CCW{binID});
        RPBDataRAllCCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CCW{binID},RPB.gainOne.dataR{flyID}.CCW{binID}),RPB.gainTwo.dataR{flyID}.CCW{binID});
        allPBDataRCCW = vertcat(LPBDataRAllCCW, RPBDataRAllCCW);
        LPBDataGAllCCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CCW{binID},LPB.gainOne.dataG{flyID}.CCW{binID}),LPB.gainTwo.dataG{flyID}.CCW{binID});
        RPBDataGAllCCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CCW{binID},RPB.gainOne.dataG{flyID}.CCW{binID}),RPB.gainTwo.dataG{flyID}.CCW{binID});
        allPBDataGCCW = vertcat(LPBDataGAllCCW, RPBDataGAllCCW);
        
        if ~isempty(LPBDataRAllCCW)
            LRCCW = mean(LPBDataRAllCCW,2);
            LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
            RRCCW = mean(RPBDataRAllCCW,2);
            RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
            RCCW = vertcat(LRCCW,RRCCW);
            LGCCW = mean(LPBDataGAllCCW,2);
            LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
            RGCCW = mean(RPBDataGAllCCW,2);
            RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
            GCCW = vertcat(LGCCW,RGCCW);

            LRPk = find(LRCCW == max(LRCCW));
            LGPk = find(LGCCW == max(LGCCW));
            pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
            RRPk = find(RRCCW == max(RRCCW));
            RGPk = find(RGCCW == max(RGCCW));
            pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
        
            subplot(4*numFlies,plotRange+1,1+binID+4*(plotRange+1)*(flyID-1)+2*(plotRange+1))
            LBin = mean(LPBDataRAllCCW,2);
            RBin = mean(RPBDataRAllCCW,2);
            imagesc(vertcat(LBin,RBin)'-1);
            line([(redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
                (redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            line([(11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
                (11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            text(2,0,num2str(pkDiffL));
            text(11,0,num2str(pkDiffR));
            axis off;
            caxis([0 0.5]);

            subplot(4*numFlies,plotRange+1,1+binID+4*(plotRange+1)*(flyID-1)+3*(plotRange+1))
            LBin = mean(LPBDataGAllCCW,2);
            RBin = mean(RPBDataGAllCCW,2);
            imagesc(vertcat(LBin,RBin)'-1);
            line([(greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
                (greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            line([(11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
                (11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
                [0.5 1.5],'Color','k','LineWidth',2);
            axis off;
            caxis([0 0.5]);
        end
    end
    
end

colormap(brewermap(64, 'Blues'));
blues = colormap(brewermap(64, 'Blues'));
greens(:,1) = blues(:,1);
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
colormap(greens);
% magentas = blues;
% magentas(:,1) = blues(:,3);
% magentas(:,2) = blues(:,1);
% magentas(:,3) = blues(:,3);
% colormap(magentas);

set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PBFig,strcat('All_vRotVsPBPosAll_OnePeak_NormColor_RAlign_greens'),'-dpdf');

%% Plot the above for a given fly, pulling out the two directions separately, using cross sections.
flyID = 4;

bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);

num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

plotRange = 180/vRSpan;
for binID=1:plotRange
    LPBDataRAllCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CW{binID},LPB.gainOne.dataR{flyID}.CW{binID}),LPB.gainTwo.dataR{flyID}.CW{binID});
    RPBDataRAllCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CW{binID},RPB.gainOne.dataR{flyID}.CW{binID}),RPB.gainTwo.dataR{flyID}.CW{binID});
    LPBDataGAllCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CW{binID},LPB.gainOne.dataG{flyID}.CW{binID}),LPB.gainTwo.dataG{flyID}.CW{binID});
    RPBDataGAllCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CW{binID},RPB.gainOne.dataG{flyID}.CW{binID}),RPB.gainTwo.dataG{flyID}.CW{binID});

    if ~isempty(LPBDataRAllCW)
        LRCW = mean(LPBDataRAllCW,2);
        LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
        RRCW = mean(RPBDataRAllCW,2);
        RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
        RCW = vertcat(LRCW,RRCW);
        LGCW = mean(LPBDataGAllCW,2);
        LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
        RGCW = mean(RPBDataGAllCW,2);
        RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
        GCW = vertcat(LGCW,RGCW);

        LRPk = find(LRCW == max(LRCW));
        LGPk = find(LGCW == max(LGCW));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRCW == max(RRCW));
        RGPk = find(RGCW == max(RGCW));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

        rSpanL = redSpan(1:8);
        rSpanR = redSpan(9:16);
        gSpanL = greenSpan(1:8);
        gSpanR = greenSpan(9:16);
        
               
        subplot(2,plotRange,binID)
        hold on;
        LBin = mean(LPBDataRAllCW,2);
        LStd = std(LPBDataRAllCW,[],2);
        RBin = mean(RPBDataRAllCW,2);
        RStd = std(RPBDataRAllCW,[],2);
        plot(rSpanL,LBin(rSpanL)-1,'m','LineWidth',2);
        plot(rSpanR,RBin(gSpanL)-1,'m','LineWidth',2);
        plot(rSpanL,LBin(rSpanL)-LStd(rSpanL)-1,'m','LineWidth',1);
        plot(rSpanR,RBin(gSpanL)-RStd(gSpanL)-1,'m','LineWidth',1);
        plot(rSpanL,LBin(rSpanL)+LStd(rSpanL)-1,'m','LineWidth',1);
        plot(rSpanR,RBin(gSpanL)+RStd(gSpanL)-1,'m','LineWidth',1);
        
        line([(redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
            (redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','m','LineWidth',2);
        line([(11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
            (11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','m','LineWidth',2);

        LBin = mean(LPBDataGAllCW,2);
        LStd = std(LPBDataGAllCW,[],2);
        RBin = mean(RPBDataGAllCW,2);
        RStd = std(RPBDataGAllCW,[],2);
        gSpanL = greenSpan(1:8);
        gSpanR = greenSpan(9:16);
        rSpanL = redSpan(1:8);
        rSpanR = redSpan(9:16);
        plot(gSpanL,LBin(gSpanL)-1,'g','LineWidth',2);
        plot(gSpanR,RBin(rSpanL)-1,'g','LineWidth',2);
        plot(gSpanL,LBin(gSpanL)-LStd(gSpanL)-1,'g','LineWidth',1);
        plot(gSpanR,RBin(rSpanL)-RStd(rSpanL)-1,'g','LineWidth',1);
        plot(gSpanL,LBin(gSpanL)+LStd(gSpanL)-1,'g','LineWidth',1);
        plot(gSpanR,RBin(rSpanL)+RStd(rSpanL)-1,'g','LineWidth',1);
        
        line([(greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
            (greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','g','LineWidth',2);
        line([(11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
            (11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','g','LineWidth',2);
        ylim([0 1.5]);;
    end

    LPBDataRAllCCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CCW{binID},LPB.gainOne.dataR{flyID}.CCW{binID}),LPB.gainTwo.dataR{flyID}.CCW{binID});
    RPBDataRAllCCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CCW{binID},RPB.gainOne.dataR{flyID}.CCW{binID}),RPB.gainTwo.dataR{flyID}.CCW{binID});
    allPBDataRCCW = vertcat(LPBDataRAllCCW, RPBDataRAllCCW);
    LPBDataGAllCCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CCW{binID},LPB.gainOne.dataG{flyID}.CCW{binID}),LPB.gainTwo.dataG{flyID}.CCW{binID});
    RPBDataGAllCCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CCW{binID},RPB.gainOne.dataG{flyID}.CCW{binID}),RPB.gainTwo.dataG{flyID}.CCW{binID});
    allPBDataGCCW = vertcat(LPBDataGAllCCW, RPBDataGAllCCW);

    if ~isempty(LPBDataRAllCCW)
        LRCCW = mean(LPBDataRAllCCW,2);
        LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
        RRCCW = mean(RPBDataRAllCCW,2);
        RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
        RCCW = vertcat(LRCCW,RRCCW);
        LGCCW = mean(LPBDataGAllCCW,2);
        LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
        RGCCW = mean(RPBDataGAllCCW,2);
        RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
        GCCW = vertcat(LGCCW,RGCCW);

        LRPk = find(LRCCW == max(LRCCW));
        LGPk = find(LGCCW == max(LGCCW));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRCCW == max(RRCCW));
        RGPk = find(RGCCW == max(RGCCW));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

        subplot(2,plotRange,binID+plotRange)
        hold on;
        LBin = mean(LPBDataRAllCCW,2);
        LStd = std(LPBDataRAllCCW,[],2);
        RBin = mean(RPBDataRAllCCW,2);
        RStd = std(RPBDataRAllCCW,[],2);
        plot(rSpanL,LBin(rSpanL)-1,'m','LineWidth',2);
        plot(rSpanR,RBin(gSpanL)-1,'m','LineWidth',2);
        plot(rSpanL,LBin(rSpanL)-LStd(rSpanL)-1,'m','LineWidth',1);
        plot(rSpanR,RBin(gSpanL)-RStd(gSpanL)-1,'m','LineWidth',1);
        plot(rSpanL,LBin(rSpanL)+LStd(rSpanL)-1,'m','LineWidth',1);
        plot(rSpanR,RBin(gSpanL)+RStd(gSpanL)-1,'m','LineWidth',1);
        
        line([(redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
            (redSpan(1)-1+circ_mean(angsraw,LBin(redSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','m','LineWidth',2);
        line([(11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
            (11-redSpan(1)+circ_mean(angsraw,RBin(sort(10-redSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','m','LineWidth',2);
        text(2,0,num2str(pkDiffL));
        text(11,0,num2str(pkDiffR));

        LBin = mean(LPBDataGAllCCW,2);
        LStd = std(LPBDataGAllCCW,[],2);
        RBin = mean(RPBDataGAllCCW,2);
        RStd = std(RPBDataGAllCCW,[],2);
        plot(gSpanL,LBin(gSpanL)-1,'g','LineWidth',2);
        plot(gSpanR,RBin(rSpanL)-1,'g','LineWidth',2);
        plot(gSpanL,LBin(gSpanL)-LStd(gSpanL)-1,'g','LineWidth',1);
        plot(gSpanR,RBin(rSpanL)-RStd(rSpanL)-1,'g','LineWidth',1);
        plot(gSpanL,LBin(gSpanL)+LStd(gSpanL)-1,'g','LineWidth',1);
        plot(gSpanR,RBin(rSpanL)+RStd(rSpanL)-1,'g','LineWidth',1);
        line([(greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)...
            (greenSpan(1)-1+circ_mean(angsraw,LBin(greenSpan(1:8)))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','g','LineWidth',2);
        line([(11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)...
            (11-greenSpan(1)+circ_mean(angsraw,RBin(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)+num_ROIs/2)],...
            [0 1.5],'Color','g','LineWidth',2);
        ylim([0 1.5]);
    end
end
    
set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompR,strcat('BumpCompGAlignNormPlots_Red',redLine,'_Green',greenLine,'AccIncEx_Profile'),'-dpdf');

%% Look at the PVA offset on either side across flies
glomShift = 5;
vRMin = 0;
vRMax = 720;
vRSpan = 30;

[LPB.dark.dataR, RPB.dark.dataR, LPB.dark.dataG, RPB.dark.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'dark',glomShift);
% [LPB.gainOne.dataR, RPB.gainOne.dataR, LPB.gainOne.dataG, RPB.gainOne.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainOne',glomShift);
% [LPB.gainTwo.dataR, RPB.gainTwo.dataR, LPB.gainTwo.dataG, RPB.gainTwo.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainTwo',glomShift);

num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

plotRange = 150/vRSpan;

% Calculate the PVA differences
LPVAdiff = zeros(numFlies,plotRange*2+1);
RPVAdiff = zeros(numFlies,plotRange*2+1);
LGPVAstren = zeros(numFlies,plotRange*2+1);
LRPVAstren = zeros(numFlies,plotRange*2+1);
RGPVAstren = zeros(numFlies,plotRange*2+1);
RRPVAstren = zeros(numFlies,plotRange*2+1);
for flyID = 1:length(allFlyData)
    
    LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
    RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
    LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
    RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;
            
    LBinG = mean(LPBDataGAllStop,2);
    RBinG = mean(RPBDataGAllStop,2);
    LBinR = mean(LPBDataRAllStop,2);
    RBinR = mean(RPBDataRAllStop,2);
    LPVAdiff(flyID,plotRange+1) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
    LPVAdiff(flyID,plotRange+1) = min(LPVAdiff(flyID,plotRange+1),9-LPVAdiff(flyID,plotRange+1));
    LGPVAstren(flyID,plotRange+1) = circ_r(angsraw,LBinG(greenSpan(1:8)));
    LRPVAstren(flyID,plotRange+1) = circ_r(angsraw,LBinR(redSpan(1:8)));
    RPVAdiff(flyID,plotRange+1) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
    RPVAdiff(flyID,plotRange+1) = min(RPVAdiff(flyID,plotRange+1),9-RPVAdiff(flyID,plotRange+1));
    RGPVAstren(flyID,plotRange+1) = circ_r(angsraw,RBinG(sort(10-greenSpan(1:8))));
    RRPVAstren(flyID,plotRange+1) = circ_r(angsraw,RBinR(sort(10-redSpan(1:8))));
    
    for binID=1:plotRange
        LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
        RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
        LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
        RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};
        
        if ~isempty(LPBDataRAllCW)
            LBinG = mean(LPBDataGAllCW,2);
            RBinG = mean(RPBDataGAllCW,2);
            LBinR = mean(LPBDataRAllCW,2);
            RBinR = mean(RPBDataRAllCW,2);
            LPVAdiff(flyID,plotRange+1-binID) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
            %LPVAdiff(flyID,plotRange+1-binID) = min(LPVAdiff(flyID,plotRange+1-binID),9-LPVAdiff(flyID,plotRange+1-binID));
            LGPVAstren(flyID,plotRange+1-binID) = circ_r(angsraw,LBinG(greenSpan(1:8)));
            LRPVAstren(flyID,plotRange+1-binID) = circ_r(angsraw,LBinR(redSpan(1:8)));
            RPVAdiff(flyID,plotRange+1-binID) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
            %RPVAdiff(flyID,plotRange+1-binID) = min(RPVAdiff(flyID,plotRange+1-binID),9-RPVAdiff(flyID,plotRange+1-binID));
            RGPVAstren(flyID,plotRange+1-binID) = circ_r(angsraw,RBinG(sort(10-greenSpan(1:8))));
            RRPVAstren(flyID,plotRange+1-binID) = circ_r(angsraw,RBinR(sort(10-redSpan(1:8))));
        end
        
        LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
        RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
        LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
        RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};
        
        if ~isempty(LPBDataRAllCCW)
            LBinG = mean(LPBDataGAllCCW,2);
            RBinG = mean(RPBDataGAllCCW,2);
            LBinR = mean(LPBDataRAllCCW,2);
            RBinR = mean(RPBDataRAllCCW,2);
            LPVAdiff(flyID,plotRange+binID) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
%             LPVAdiff(flyID,plotRange+binID) = min(LPVAdiff(flyID,plotRange+binID),9-LPVAdiff(flyID,plotRange+binID));
            LGPVAstren(flyID,plotRange+1+binID) = circ_r(angsraw,LBinG(greenSpan(1:8)));
            LRPVAstren(flyID,plotRange+1+binID) = circ_r(angsraw,LBinR(redSpan(1:8)));
            RPVAdiff(flyID,plotRange+binID) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
%             RPVAdiff(flyID,plotRange+binID) = min(RPVAdiff(flyID,plotRange+binID),9-RPVAdiff(flyID,plotRange+binID));
            RGPVAstren(flyID,plotRange+1+binID) = circ_r(angsraw,RBinG(sort(10-greenSpan(1:8))));
            RRPVAstren(flyID,plotRange+1+binID) = circ_r(angsraw,RBinR(sort(10-redSpan(1:8))));
        end
    end
    
end

PBFig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LPVAdiff(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RPVAdiff(:,velPlt),'r');
end
ylim([-2 2]);
subplot(2,1,2)
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LGPVAstren(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RGPVAstren(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LRPVAstren(:,velPlt),'r');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RRPVAstren(:,velPlt),'r');
end

save(strcat('PVADiff_',greenLine,'-',redLine,'2.mat'),'LPVAdiff','RPVAdiff','LGPVAstren','LRPVAstren','RGPVAstren','RRPVAstren');

%% Plot the peak width vs. the rotational velocity

% Specify the range of rotational velocities to consider
vRMin = 0;
vRMax = 720;
vRSpan = 15;

glomShift = 5;

[LPB.dark.dataR, RPB.dark.dataR, LPB.dark.dataG, RPB.dark.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'dark',glomShift);
[LPB.gainOne.dataR, RPB.gainOne.dataR, LPB.gainOne.dataG, RPB.gainOne.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainOne',glomShift);
[LPB.gainTwo.dataR, RPB.gainTwo.dataR, LPB.gainTwo.dataG, RPB.gainTwo.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainTwo',glomShift);

plotMax = 180;
plotRange = plotMax/vRSpan;

LGpkMax = zeros(length(allFlyData),1+2*plotRange);
LGpkWid = zeros(length(allFlyData),1+2*plotRange);
RGpkMax = zeros(length(allFlyData),1+2*plotRange);
RGpkWid = zeros(length(allFlyData),1+2*plotRange);
LRpkMax = zeros(length(allFlyData),1+2*plotRange);
LRpkWid = zeros(length(allFlyData),1+2*plotRange);
RRpkMax = zeros(length(allFlyData),1+2*plotRange);
RRpkWid = zeros(length(allFlyData),1+2*plotRange);

% Look at the red channel   
for flyID = 1:length(allFlyData)
    
    LPBDataRAllStop = horzcat(horzcat(LPB.dark.dataR{flyID}.Stop,LPB.gainOne.dataR{flyID}.Stop),LPB.gainTwo.dataR{flyID}.Stop);
    RPBDataRAllStop = horzcat(horzcat(RPB.dark.dataR{flyID}.Stop,RPB.gainOne.dataR{flyID}.Stop),RPB.gainTwo.dataR{flyID}.Stop);
    LPBDataGAllStop = horzcat(horzcat(LPB.dark.dataG{flyID}.Stop,LPB.gainOne.dataG{flyID}.Stop),LPB.gainTwo.dataG{flyID}.Stop);
    RPBDataGAllStop = horzcat(horzcat(RPB.dark.dataG{flyID}.Stop,RPB.gainOne.dataG{flyID}.Stop),RPB.gainTwo.dataG{flyID}.Stop);
    
    LGStop = mean(LPBDataGAllStop(greenSpan(1:8),:),2);
    LGpkMax(flyID, 1+plotRange) = max(LGStop)-min(LGStop);
    LGpkRng = find(LGStop-min(LGStop)>0.5*LGpkMax(flyID, 1+plotRange));
    if ~isempty(LGpkRng)
       LGpkWid(flyID, 1+plotRange) = LGpkRng(end)-LGpkRng(1);
    end

    RGStop = mean(RPBDataGAllStop(redSpan(1:8),:),2);
    RGpkMax(flyID, 1+plotRange) = max(RGStop)-min(RGStop);
    RGpkRng = find(RGStop-min(RGStop)>0.5*RGpkMax(flyID, 1+plotRange));
    if ~isempty(RGpkRng)
       RGpkWid(flyID, 1+plotRange) = RGpkRng(end)-RGpkRng(1);
    end
    
    LRStop = mean(LPBDataRAllStop(redSpan(1:8),:),2);
    LRpkMax(flyID, 1+plotRange) = max(LRStop)-min(LRStop);
    LRpkRng = find(LRStop-min(LRStop)>0.5*LRpkMax(flyID, 1+plotRange));
    if ~isempty(LRpkRng)
       LRpkWid(flyID, 1+plotRange) = LRpkRng(end)-LRpkRng(1);
    end

    RRStop = mean(RPBDataRAllStop(greenSpan(1:8),:),2);
    RRpkMax(flyID, 1+plotRange) = max(RRStop)-min(RRStop);
    RRpkRng = find(RRStop-min(RRStop)>0.5*RRpkMax(flyID, 1+plotRange));
    if ~isempty(RRpkRng)
       RRpkWid(flyID, 1+plotRange) = RRpkRng(end)-RRpkRng(1);
    end
    
    for binID=1:plotRange
        LPBDataRAllCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CW{binID},LPB.gainOne.dataR{flyID}.CW{binID}),LPB.gainTwo.dataR{flyID}.CW{binID});
        RPBDataRAllCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CW{binID},RPB.gainOne.dataR{flyID}.CW{binID}),RPB.gainTwo.dataR{flyID}.CW{binID});
        LPBDataGAllCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CW{binID},LPB.gainOne.dataG{flyID}.CW{binID}),LPB.gainTwo.dataG{flyID}.CW{binID});
        RPBDataGAllCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CW{binID},RPB.gainOne.dataG{flyID}.CW{binID}),RPB.gainTwo.dataG{flyID}.CW{binID});

        if ~isempty(LPBDataRAllCW)
            LGStop = mean(LPBDataGAllCW(greenSpan(1:8),:),2);
            LGpkMax(flyID, 1+plotRange-binID) = max(LGStop)-min(LGStop);
            LGpkRng = find(LGStop-min(LGStop)>0.5*LGpkMax(flyID, 1+plotRange-binID));
            if ~isempty(LGpkRng)
               LGpkWid(flyID, 1+plotRange-binID) = LGpkRng(end)-LGpkRng(1);
            end

            RGStop = mean(RPBDataGAllCW(redSpan(1:8),:),2);
            RGpkMax(flyID, 1+plotRange-binID) = max(RGStop)-min(RGStop);
            RGpkRng = find(RGStop-min(RGStop)>0.5*RGpkMax(flyID, 1+plotRange-binID));
            if ~isempty(RGpkRng)
               RGpkWid(flyID, 1+plotRange-binID) = RGpkRng(end)-RGpkRng(1);
            end

            LRStop = mean(LPBDataRAllCW(redSpan(1:8),:),2);
            LRpkMax(flyID, 1+plotRange-binID) = max(LRStop)-min(LRStop);
            LRpkRng = find(LRStop-min(LRStop)>0.5*LRpkMax(flyID, 1+plotRange-binID));
            if ~isempty(LRpkRng)
               LRpkWid(flyID, 1+plotRange-binID) = LRpkRng(end)-LRpkRng(1);
            end

            RRStop = mean(RPBDataRAllCW(greenSpan(1:8),:),2);
            RRpkMax(flyID, 1+plotRange-binID) = max(RRStop)-min(RRStop);
            RRpkRng = find(RRStop-min(RRStop)>0.5*RRpkMax(flyID, 1+plotRange-binID));
            if ~isempty(RRpkRng)
               RRpkWid(flyID, 1+plotRange-binID) = RRpkRng(end)-RRpkRng(1);
            end
        end
        
        LPBDataRAllCCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CCW{binID},LPB.gainOne.dataR{flyID}.CCW{binID}),LPB.gainTwo.dataR{flyID}.CCW{binID});
        RPBDataRAllCCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CCW{binID},RPB.gainOne.dataR{flyID}.CCW{binID}),RPB.gainTwo.dataR{flyID}.CCW{binID});
        LPBDataGAllCCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CCW{binID},LPB.gainOne.dataG{flyID}.CCW{binID}),LPB.gainTwo.dataG{flyID}.CCW{binID});
        RPBDataGAllCCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CCW{binID},RPB.gainOne.dataG{flyID}.CCW{binID}),RPB.gainTwo.dataG{flyID}.CCW{binID});

        if ~isempty(LPBDataRAllCCW)
            LGStop = mean(LPBDataGAllCCW(greenSpan(1:8),:),2);
            LGpkMax(flyID, 1+plotRange+binID) = max(LGStop)-min(LGStop);
            LGpkRng = find(LGStop-min(LGStop)>0.5*LGpkMax(flyID, 1+plotRange+binID));
            if ~isempty(LGpkRng)
               LGpkWid(flyID, 1+plotRange+binID) = LGpkRng(end)-LGpkRng(1);
            end

            RGStop = mean(RPBDataGAllCCW(redSpan(1:8),:),2);
            RGpkMax(flyID, 1+plotRange+binID) = max(RGStop)-min(RGStop);
            RGpkRng = find(RGStop-min(RGStop)>0.5*RGpkMax(flyID, 1+plotRange+binID));
            if ~isempty(RGpkRng)
               RGpkWid(flyID, 1+plotRange+binID) = RGpkRng(end)-RGpkRng(1);
            end

            LRStop = mean(LPBDataRAllCCW(redSpan(1:8),:),2);
            LRpkMax(flyID, 1+plotRange+binID) = max(LRStop)-min(LRStop);
            LRpkRng = find(LRStop-min(LRStop)>0.5*LRpkMax(flyID, 1+plotRange+binID));
            if ~isempty(LRpkRng)
               LRpkWid(flyID, 1+plotRange+binID) = LRpkRng(end)-LRpkRng(1);
            end

            RRStop = mean(RPBDataRAllCCW(greenSpan(1:8),:),2);
            RRpkMax(flyID, 1+plotRange+binID) = max(RRStop)-min(RRStop);
            RRpkRng = find(RRStop-min(RRStop)>0.5*RRpkMax(flyID, 1+plotRange+binID));
            if ~isempty(RRpkRng)
               RRpkWid(flyID, 1+plotRange+binID) = RRpkRng(end)-RRpkRng(1);
            end
        end
    end
end

pkWid = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID=1:length(allFlyData)
    subplot(2,2,1);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,LGpkWid(flyID,:));
    
    subplot(2,2,2);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,RGpkWid(flyID,:));
    
    subplot(2,2,3);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,LRpkWid(flyID,:));
    
    subplot(2,2,4);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,RRpkWid(flyID,:));
end

pkMax = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID=1:length(allFlyData)
    subplot(2,2,1);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,LGpkMax(flyID,:));
    
    subplot(2,2,2);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,RGpkMax(flyID,:));
    
    subplot(2,2,3);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,LRpkMax(flyID,:));
    
    subplot(2,2,4);
    hold on;
    scatter(-plotMax:vRSpan:plotMax,RRpkMax(flyID,:));
end


%% Align the peaks to one another as sorted by side and velocity
vRMin = 0;
vRMax = 360;
vRSpan = 30;
[LPB.dark.dataR, LPB.dark.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'dark');
[LPB.gainOne.dataR, LPB.gainOne.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'gainOne');
[LPB.gainTwo.dataR, LPB.gainTwo.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'gainTwo');
[RPB.dark.dataR, RPB.dark.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'dark');
[RPB.gainOne.dataR, RPB.gainOne.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'gainOne');
[RPB.gainTwo.dataR, RPB.gainTwo.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'gainTwo');

vRBinNum = round((vRMax-vRMin)/vRSpan);

for flyID = 1:length(allFlyData)
    PBFig = figure('units','normalized','outerposition',[0 0 1 1]);
    % Plot the dark data
    subplot(6,vRBinNum+1,1)
    hold on;
    if ~isempty(LPB.dark.dataR{flyID}.Stop)
        plot(1:8,circshift(median(LPB.dark.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(2:9,circshift(median(LPB.dark.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    title(strcat('Left PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(LPB.dark.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,binID+1);
            hold on;
            plot(1:8,circshift(median(LPB.dark.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.dark.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(LPB.dark.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,binID+vRBinNum+2);
            hold on;
            plot(1:8,circshift(median(LPB.dark.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.dark.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.dark.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.dark.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    
    % Plot the gainOne data
    subplot(6,vRBinNum+1,2*(vRBinNum+1)+1)
    hold on;
    if ~isempty(LPB.gainOne.dataR{flyID}.Stop)
        plot(1:8,circshift(median(LPB.gainOne.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(2:9,circshift(median(LPB.gainOne.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    ylabel('gainOne');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(LPB.gainOne.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,2*(vRBinNum+1)+binID+1);
            hold on;
            plot(1:8,circshift(median(LPB.gainOne.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.gainOne.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(LPB.gainOne.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,3*(vRBinNum+1)+binID+1);
            hold on;
            plot(1:8,circshift(median(LPB.gainOne.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.gainOne.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.gainOne.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.gainOne.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    
    % Plot the gainTwo data
    subplot(6,vRBinNum+1,4*(vRBinNum+1)+1)
    hold on;
    if ~isempty(LPB.gainTwo.dataR{flyID}.Stop)
        plot(1:8,circshift(median(LPB.gainTwo.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(2:9,circshift(median(LPB.gainTwo.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    ylabel('gainTwo');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(LPB.gainTwo.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,4*(vRBinNum+1)+binID+1);
            hold on;
            plot(1:8,circshift(median(LPB.gainTwo.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.gainTwo.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(LPB.gainTwo.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,5*(vRBinNum+1)+binID+1);
            hold on;
            plot(1:8,circshift(median(LPB.gainTwo.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(1:8,circshift(quantile(LPB.gainTwo.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(2:9,circshift(median(LPB.gainTwo.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(2:9,circshift(quantile(LPB.gainTwo.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    
    set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPos_Left'),'-dpdf');
end

for flyID = 1:length(allFlyData)
    PBFig = figure('units','normalized','outerposition',[0 0 1 1]);
    % Plot the dark data
    subplot(6,vRBinNum+1,1)
    hold on;
    if ~isempty(RPB.dark.dataR{flyID}.Stop)
        plot(11:18,circshift(median(RPB.dark.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(10:17,circshift(median(RPB.dark.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    title(strcat('Right PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(RPB.dark.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,binID+1);
            hold on;
            plot(11:18,circshift(median(RPB.dark.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.dark.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(RPB.dark.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,binID+vRBinNum+2);
            hold on;
            plot(11:18,circshift(median(RPB.dark.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.dark.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.dark.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.dark.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    
    % Plot the gainOne data
    subplot(6,vRBinNum+1,2*(vRBinNum+1)+1)
    hold on;
    if ~isempty(RPB.gainOne.dataR{flyID}.Stop)
        plot(11:18,circshift(median(RPB.gainOne.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(10:17,circshift(median(RPB.gainOne.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    ylabel('gainOne');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(RPB.gainOne.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,2*(vRBinNum+1)+binID+1);
            hold on;
            plot(11:18,circshift(median(RPB.gainOne.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.gainOne.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(RPB.gainOne.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,3*(vRBinNum+1)+binID+1);
            hold on;
            plot(11:18,circshift(median(RPB.gainOne.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.gainOne.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.gainOne.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.gainOne.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    
    % Plot the gainTwo data
    subplot(6,vRBinNum+1,4*(vRBinNum+1)+1)
    hold on;
    if ~isempty(RPB.gainTwo.dataR{flyID}.Stop)
        plot(11:18,circshift(median(RPB.gainTwo.dataR{flyID}.Stop,2),4),'r','LineWidth',2);
        plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.Stop,0.25,2),4),'r','LineWidth',1);
        plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.Stop,0.75,2),4),'r','LineWidth',1);
        plot(10:17,circshift(median(RPB.gainTwo.dataG{flyID}.Stop,2),4),'g','LineWidth',2);
        plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.Stop,0.25,2),4),'g','LineWidth',1);
        plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.Stop,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    ylabel('gainTwo');
    xlabel('0');
    for binID=1:vRBinNum
        if ~isempty(RPB.gainTwo.dataR{flyID}.CW{binID})
            subplot(6,vRBinNum+1,4*(vRBinNum+1)+binID+1);
            hold on;
            plot(11:18,circshift(median(RPB.gainTwo.dataR{flyID}.CW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.CW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.CW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.gainTwo.dataG{flyID}.CW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.CW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.CW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
        if ~isempty(RPB.gainTwo.dataR{flyID}.CCW{binID})
            subplot(6,vRBinNum+1,5*(vRBinNum+1)+binID+1);
            hold on;
            plot(11:18,circshift(median(RPB.gainTwo.dataR{flyID}.CCW{binID},2),4),'r','LineWidth',2);
            plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.CCW{binID},0.25,2),4),'r','LineWidth',1);
            plot(11:18,circshift(quantile(RPB.gainTwo.dataR{flyID}.CCW{binID},0.75,2),4),'r','LineWidth',1);
            plot(10:17,circshift(median(RPB.gainTwo.dataG{flyID}.CCW{binID},2),4),'g','LineWidth',2);
            plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.CCW{binID},0.25,2),4),'g','LineWidth',1);
            plot(10:17,circshift(quantile(RPB.gainTwo.dataG{flyID}.CCW{binID},0.75,2),4),'g','LineWidth',1);
            ylim([1 2]);
        end
    end
    set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPos_Right'),'-dpdf');
end

%% Plot a summary for individual dark, 1x, and 2x trials - amplitude and vrot

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID});
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

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
            vRot = smooth(diff(OffsetRotMatchUnwrap),11);

            % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            hold on;
            plot(OffsetRotMatch(:,1),max(RROIaveMaxFilt(1:9,minFG:maxFG),[],1),'r');
            plot(OffsetRotMatch(:,1),max(RROIaveMaxFilt(10:18,minFG:maxFG),[],1),'--r');
            title(redLine);

            subplot(4,1,2);
            hold on;
            plot(OffsetRotMatch(:,1),max(GROIaveMaxFilt(1:9,minFG:maxFG),[],1),'g');
            plot(OffsetRotMatch(:,1),max(GROIaveMaxFilt(10:18,minFG:maxFG),[],1),'--g');
            title(greenLine);

            subplot(4,1,3);
            plot(OffsetRotMatch(1:end-1,1),vRot,'k');
            ylabel('rotational speed (rad/sec)');
            
            subplot(4,1,4);
            hold on;
            Rsig = smooth(max(RROIaveMaxFilt(:,minFG:maxFG),[],1),11);
            plot(OffsetRotMatch(:,1),(Rsig-min(Rsig))./(max(Rsig-min(Rsig))),'r');
            Gsig = smooth(max(GROIaveMaxFilt(:,minFG:maxFG),[],1),11);
            plot(OffsetRotMatch(:,1),(Gsig-min(Gsig))./(max(Gsig-min(Gsig))),'g');
            plot(OffsetRotMatch(1:end-1,1),abs(vRot)./max(abs(vRot)),'k');
            ylim([0 1]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_Dark_Summary_',num2str(darkID)),'-dpdf');
            print(RGcomp,strcat('Fly',num2str(flyID),'_Dark_Summary_AmpVrot_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);
            
            
             % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);
            
            RVCorr = zeros(201,1);
            GVCorr = zeros(201,1);
            RGCorr = zeros(201,1);
            speed = abs(vRot);
                
            for i=1:floor(length(RVCorr)/2)
                Rccs = corrcoef(Rsig(i:end-1),speed(1:end-i+1));
                Gccs = corrcoef(Gsig(i:end-1),speed(1:end-i+1));
                RGccs = corrcoef(Gsig(i:end),Rsig(1:end-i+1));
                RVCorr(floor(length(RVCorr)/2)+i-1) = Rccs(2,1);
                GVCorr(floor(length(RVCorr)/2)+i-1) = Gccs(2,1);
                RGCorr(floor(length(RGCorr)/2)+i-1) = RGccs(2,1);
            end
            
            for i=1:floor(length(RVCorr)/2)
                Rccs = corrcoef(Rsig(1:end-i),speed(i:end));
                Gccs = corrcoef(Gsig(1:end-i),speed(i:end));
                RGccs = corrcoef(Gsig(1:end-i+1),Rsig(i:end));
                RVCorr(floor(length(RVCorr)/2)-i+1) = Rccs(2,1);
                GVCorr(floor(length(RVCorr)/2)-i+1) = Gccs(2,1);
                RGCorr(floor(length(RGCorr)/2)-i+1) = RGccs(2,1);
            end
            
            tStep = mean(diff(OffsetRotMatch(:,1)));

            subplot(4,1,1);
            plot([-100:100]*tStep,fliplr(RVCorr),'r');
            title(redLine);
            ylim([0 1]);
            xlim([-100*tStep 100*tStep]);
            
            subplot(4,1,2);
            plot([-100:100]*tStep,fliplr(GVCorr),'g');
            title(greenLine);
            ylim([0 1]);
            xlim([-100*tStep 100*tStep]);
            
            subplot(4,1,3);
            hold on;
            plot([-100:100]*tStep,fliplr(RVCorr),'r');
            plot([-100:100]*tStep,fliplr(GVCorr),'g');
            ylim([0 1]);
            xlim([-100*tStep 100*tStep]);
            
            subplot(4,1,4);
            plot([-100:100]*tStep,fliplr(RGCorr),'k');
            title(strcat(redLine,' vs ',greenLine));
            ylim([0 1]);
            xlim([-100*tStep 100*tStep]);
            
            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
            print(RGcomp,strcat('Fly',num2str(flyID),'_Dark_Summary_AmpVrotCorr_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);

            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
   
end

%% Create array of bump shapes as sorted by velocity and sign of acceleration (for dark case)
% number of frames over which to S-G filter
framesSmooth = 11;
tDelay = 0;

% ROIrange = [1:9]; side = 'left';
ROIrange = [10:18]; side = 'right';
shift = 5;

% Threshold velocities to consider
% velThreshLC = [0:5:90];
% velThreshUC = [5:5:95];
velThreshLC = [0:20:80];
velThreshUC = [20:20:100];
XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID})
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = framesSmooth;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

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
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
%             OffsetRotUnwrap = sgolayfilt(OffsetRotUnwrap,sgolayOrder,sgolayWindow);
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,framesSmooth);
                    
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins by increasing or
                % decreasing velocity
                RturnsGinc = [];
                RturnsGdec = [];
                RturnsRinc = [];
                RturnsRdec = [];
                LturnsGinc = [];
                LturnsGdec = [];
                LturnsRinc = [];
                LturnsRdec = [];

                velThreshUCNow = velThreshUC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                velThreshLCNow = velThreshLC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                aFly = diff(diff(OffsetRotUnwrap));

                for i=minFG+1+tDelay:maxFG-2
                    if aFly(i+1-minFG-tDelay) > 0 
                        if OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) > velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) < velThreshUCNow
                            RturnsGinc = horzcat(RturnsGinc, GROIaveMaxFilt(:,i));
                            RturnsRinc = horzcat(RturnsRinc, RROIaveMaxFilt(:,i));
                        elseif OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) < -velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) > -velThreshUCNow
                            LturnsGinc = horzcat(LturnsGinc, GROIaveMaxFilt(:,i));
                            LturnsRinc = horzcat(LturnsRinc, RROIaveMaxFilt(:,i));
                        end
                    elseif aFly(i+1-minFG-tDelay) < 0
                        if OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) > velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) < velThreshUCNow
                            RturnsGdec = horzcat(RturnsGdec, GROIaveMaxFilt(:,i));
                            RturnsRdec = horzcat(RturnsRdec, RROIaveMaxFilt(:,i));
                        elseif OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) < -velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG-tDelay) - OffsetRotUnwrap(i-minFG-tDelay) > -velThreshUCNow
                            LturnsGdec = horzcat(LturnsGdec, GROIaveMaxFilt(:,i));
                            LturnsRdec = horzcat(LturnsRdec, RROIaveMaxFilt(:,i));
                        end
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsGinc = RturnsGinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsGdec = RturnsGdec;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsRinc = RturnsRinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsRdec = RturnsRdec;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsGinc = LturnsGinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsGdec = LturnsGdec;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsRinc = LturnsRinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsRdec = LturnsRdec;
            end
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

% Pool together the turns across trials for ROIs 1:9
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       RturnsGallGalignInc = [];
       RturnsGallGalignDec = [];
       RturnsRallGalignInc = [];
       RturnsRallGalignDec = [];
       LturnsGallGalignInc = [];
       LturnsGallGalignDec = [];
       LturnsRallGalignInc = [];
       LturnsRallGalignDec = [];
       
       RturnsGallRalignInc = [];
       RturnsGallRalignDec = [];
       RturnsRallRalignInc = [];
       RturnsRallRalignDec = [];
       LturnsGallRalignInc = [];
       LturnsGallRalignDec = [];
       LturnsRallRalignInc = [];
       LturnsRallRalignDec = [];
        for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
            RturnsGGalignInc = [];
            RturnsGGalignDec = [];
            RturnsRGalignInc = [];
            RturnsRGalignDec = [];
            LturnsGGalignInc = [];
            LturnsGGalignDec = [];
            LturnsRGalignInc = [];
            LturnsRGalignDec = [];
            
            RturnsGRalignInc = [];
            RturnsGRalignDec = [];
            RturnsRRalignInc = [];
            RturnsRRalignDec = [];
            LturnsGRalignInc = [];
            LturnsGRalignDec = [];
            LturnsRRalignInc = [];
            LturnsRRalignDec = [];
            
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.dark{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i)));
%                     pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i)));
%                     pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i)));
                    RturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i),shift-pkPosGalign);
                    RturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i),shift-pkPosGalign);
                    RturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(ROIrange,i),shift-pkPosRalign);
                    RturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(ROIrange,i),shift-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i)));
%                     pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i)));
%                     pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i)));
                    LturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i),shift-pkPosGalign);
                    LturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i),shift-pkPosGalign);
                    LturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(ROIrange,i),shift-pkPosRalign);
                    LturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(ROIrange,i),shift-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i)));
%                     pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i)));
%                     pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i)));
                    RturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i),shift-pkPosGalign);
                    RturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i),shift-pkPosGalign);
                    RturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(ROIrange,i),shift-pkPosRalign);
                    RturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(ROIrange,i),shift-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i)));
%                     pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i)));
%                     pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i)==...
%                         min(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i)));
                    LturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i),shift-pkPosGalign);
                    LturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i),shift-pkPosGalign);
                    LturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(ROIrange,i),shift-pkPosRalign);
                    LturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(ROIrange,i),shift-pkPosRalign);
                end
            end
            
            RturnsGallGalignInc = horzcat(RturnsGallGalignInc,RturnsGGalignInc);
            RturnsRallGalignInc = horzcat(RturnsRallGalignInc,RturnsRGalignInc);
            LturnsGallGalignInc = horzcat(LturnsGallGalignInc,LturnsGGalignInc);
            LturnsRallGalignInc = horzcat(LturnsRallGalignInc,LturnsRGalignInc);
            RturnsGallGalignDec = horzcat(RturnsGallGalignDec,RturnsGGalignDec);
            RturnsRallGalignDec = horzcat(RturnsRallGalignDec,RturnsRGalignDec);
            LturnsGallGalignDec = horzcat(LturnsGallGalignDec,LturnsGGalignDec);
            LturnsRallGalignDec = horzcat(LturnsRallGalignDec,LturnsRGalignDec);
            
            RturnsGallRalignInc = horzcat(RturnsGallRalignInc,RturnsGRalignInc);
            RturnsRallRalignInc = horzcat(RturnsRallRalignInc,RturnsRRalignInc);
            LturnsGallRalignInc = horzcat(LturnsGallRalignInc,LturnsGRalignInc);
            LturnsRallRalignInc = horzcat(LturnsRallRalignInc,LturnsRRalignInc);
            RturnsGallRalignDec = horzcat(RturnsGallRalignDec,RturnsGRalignDec);
            RturnsRallRalignDec = horzcat(RturnsRallRalignDec,RturnsRRalignDec);
            LturnsGallRalignDec = horzcat(LturnsGallRalignDec,LturnsGRalignDec);
            LturnsRallRalignDec = horzcat(LturnsRallRalignDec,LturnsRRalignDec);

       end
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc = RturnsGallGalignInc;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc = RturnsRallGalignInc;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc = LturnsGallGalignInc;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc = LturnsRallGalignInc;
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec = RturnsGallGalignDec;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec = RturnsRallGalignDec;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec = LturnsGallGalignDec;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec = LturnsRallGalignDec;
       
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc = RturnsGallRalignInc;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc = RturnsRallRalignInc;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc = LturnsGallRalignInc;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc = LturnsRallRalignInc;
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec = RturnsGallRalignDec;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec = RturnsRallRalignDec;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec = LturnsGallRalignDec;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec = LturnsRallRalignDec;
    end 
end

bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec);
           RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec);
           LturnsG = horzcat(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc,XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec);
           LturnsR = horzcat(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc,XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec);
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(8,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(8,0.8,num2str(size(LnetR,2)),'Color','m');
           plot(median(LnetG,2),'--g','LineWidth',1);
           
           xlim([1 9]);
           ylim([0 2]);
           set(gca,'FontSize',14);
       end
       if flyID == 1
           title(strcat('v_{thresh} = ',num2str(XCVals{flyID}.vel{velCut}.threshold),'deg/sec'));
       end
       if flyID == numFlies
           xlabel('ROI');
       end
    end
    subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+1);
    ylabel(allFlyData{flyID}.ID);
end


set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompR,strcat(side,'PBLeadLag_Galign_tD',num2str(tDelay)),'-dpdf');

bumpCompG = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec);
           RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec);
           LturnsG = horzcat(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec);
           LturnsR = horzcat(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec);
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetG,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetG,0.75,2),'b','LineWidth',1);
           plot(median(RnetG,2),'b','LineWidth',4);
           text(8,1.1,num2str(size(RnetG,2)),'Color','b');
           plot(median(RnetR,2),'--r','LineWidth',1);
           
           plot(quantile(LnetG,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetG,0.75,2),'m','LineWidth',1);
           plot(median(LnetG,2),'m','LineWidth',4);
           text(8,0.8,num2str(size(LnetG,2)),'Color','m');
           plot(median(LnetR,2),'--r','LineWidth',1);
           
           xlim([1 9]);
           ylim([0 2]);
           set(gca,'FontSize',14);
       end
       if flyID == 1
           title(strcat('v_{thresh} = ',num2str(XCVals{flyID}.vel{velCut}.threshold),'deg/sec'));
       end
       if flyID == numFlies
           xlabel('ROI');
       end
    end
    subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+1);
    ylabel(allFlyData{flyID}.ID);
end


set(bumpCompG,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompG,strcat(side,'PBLeadLag_Ralign_tD',num2str(tDelay)),'-dpdf');

%% Plot the above for a given fly, pulling out the two directions separately.
% flyID = 7;
% flyID = 5;
flyID = 1;

bumpCompRinc = figure('units','normalized','outerposition',[0 0 1 1]);

for velCut = 1:length(velThreshLC)       
   for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;
       subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
       hold on;

       subplot(2,length(velThreshLC),velCut);
       hold on;

       plot(quantile(RnetR,0.25,2),'r','LineWidth',1);
       plot(quantile(RnetR,0.75,2),'r','LineWidth',1);
       plot(median(RnetR,2),'r','LineWidth',4);
       plot(quantile(RnetG,0.25,2),'g','LineWidth',1);
       plot(quantile(RnetG,0.75,2),'g','LineWidth',1);
       plot(median(RnetG,2),'g','LineWidth',4);

       xlim([1 16]);
       ylim([-0.5 2]);
       set(gca,'FontSize',14);


       subplot(2,length(velThreshLC),velCut+length(velThreshLC));
       hold on;

       plot(quantile(LnetR,0.25,2),'r','LineWidth',1);
       plot(quantile(LnetR,0.75,2),'r','LineWidth',1);
       plot(median(LnetR,2),'r','LineWidth',4);
       plot(quantile(LnetG,0.25,2),'g','LineWidth',1);
       plot(quantile(LnetG,0.75,2),'g','LineWidth',1);
       plot(median(LnetG,2),'g','LineWidth',4);

       xlim([1 16]);
       ylim([-0.5 2]);
       set(gca,'FontSize',14);
   end
end


bumpCompRdec = figure('units','normalized','outerposition',[0 0 1 1]);
for velCut = 1:length(velThreshLC)       
   for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;
       subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
       hold on;

       subplot(2,length(velThreshLC),velCut);
       hold on;

       plot(quantile(RnetR,0.25,2),'r','LineWidth',1);
       plot(quantile(RnetR,0.75,2),'r','LineWidth',1);
       plot(median(RnetR,2),'r','LineWidth',4);
       plot(quantile(RnetG,0.25,2),'g','LineWidth',1);
       plot(quantile(RnetG,0.75,2),'g','LineWidth',1);
       plot(median(RnetG,2),'g','LineWidth',4);

       xlim([1 16]);
       ylim([-0.5 2]);
       set(gca,'FontSize',14);


       subplot(2,length(velThreshLC),velCut+length(velThreshLC));
       hold on;

       plot(quantile(LnetR,0.25,2),'r','LineWidth',1);
       plot(quantile(LnetR,0.75,2),'r','LineWidth',1);
       plot(median(LnetR,2),'r','LineWidth',4);
       plot(quantile(LnetG,0.25,2),'g','LineWidth',1);
       plot(quantile(LnetG,0.75,2),'g','LineWidth',1);
       plot(median(LnetG,2),'g','LineWidth',4);

       xlim([1 16]);
       ylim([-0.5 2]);
       set(gca,'FontSize',14);
   end
end

set(bumpCompRinc,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompRinc,strcat(side,'PB_BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccIncEx_tD',num2str(tDelay)),'-dpdf');
set(bumpCompRdec,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompRdec,strcat(side,'PB_BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccDecEx_tD',num2str(tDelay)),'-dpdf');

%% Get the peak heights and widths 
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=1:3
    for velCut = 1:length(velThreshLC)       
       RGturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec);
       LGturnsG = horzcat(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec);
       
       RGturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec);
       LGturnsR = horzcat(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec);
       
       RRturnsG = ...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec;
       LRturnsG = ...
           XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec;
       
       RRturnsR = ...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec;
       LRturnsR = ...
           XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec;

       % Calculate the PVA differences
       RGpkMax(flyID, velCut) = max(median(RRturnsG,2))-min(median(RRturnsG,2));
       RGpkRng = find(median(RRturnsG,2)-min(median(RRturnsG,2))>0.5*RGpkMax(flyID, velCut));
       RGpkWid(flyID, velCut) = RGpkRng(end)-RGpkRng(1);
       LGpkMax(flyID, velCut) = max(median(LRturnsG,2))-min(median(LRturnsG,2));
       LGpkRng = find(median(LRturnsG,2)-min(median(LRturnsG,2))>0.5*LGpkMax(flyID, velCut));
       LGpkWid(flyID, velCut) = LGpkRng(end)-LGpkRng(1);
       
       RRpkMax(flyID, velCut) = max(median(RGturnsR,2))-min(median(RGturnsR,2));
       RRpkRng = find(median(RGturnsR,2)-min(median(RGturnsR,2))>0.5*RRpkMax(flyID, velCut));
       RRpkWid(flyID, velCut) = RRpkRng(end)-RRpkRng(1);
       LRpkMax(flyID, velCut) = max(median(LGturnsR,2))-min(median(LGturnsR,2));
       LRpkRng = find(median(LGturnsR,2)-min(median(LGturnsR,2))>0.5*LRpkMax(flyID, velCut));
       LRpkWid(flyID, velCut) = LRpkRng(end)-LRpkRng(1);
       
    end

end

PVAdiff = PVAdiff;

figure;
set(gcf,'Position',[50 50 900 900]);
subplot(2,2,1);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+20*i,squeeze(RGpkMax(:,i)),'filled','k'); 
   line([20*i-2 20*i+2],[mean(RGpkMax(:,i)) mean(RGpkMax(:,i))],'Color','k','LineWidth',2);
   scatter(zeros(numFlies,1)+20*i,squeeze(LGpkMax(:,i)),'filled','c'); 
   line([20*i-2 20*i+2],[mean(LGpkMax(:,i)) mean(LGpkMax(:,i))],'Color','c','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 95]);
xlabel('vr bins (deg/sec)');

subplot(2,2,2);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+20*i,40*squeeze(RGpkWid(:,i)),'filled','k'); 
   line([20*i-2 20*i+2],[40*mean(RGpkWid(:,i)) 40*mean(RGpkWid(:,i))],'Color','k','LineWidth',2);
   scatter(zeros(numFlies,1)+20*i,40*squeeze(LGpkWid(:,i)),'filled','c'); 
   line([20*i-2 20*i+2],[40*mean(LGpkWid(:,i)) 40*mean(LGpkWid(:,i))],'Color','c','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 95]);
xlabel('vr bins (deg/sec)');
ylim([0 250]);

subplot(2,2,3);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+20*i,squeeze(RRpkMax(:,i)),'filled','k'); 
   line([20*i-2 20*i+2],[mean(RRpkMax(:,i)) mean(RRpkMax(:,i))],'Color','k','LineWidth',2);
   scatter(zeros(numFlies,1)+20*i,squeeze(LRpkMax(:,i)),'filled','c'); 
   line([20*i-2 20*i+2],[mean(LRpkMax(:,i)) mean(LRpkMax(:,i))],'Color','c','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 95]);
xlabel('vr bins (deg/sec)');

subplot(2,2,4);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+20*i,40*squeeze(RRpkWid(:,i)),'filled','k'); 
   line([20*i-2 20*i+2],[40*mean(RRpkWid(:,i)) 40*mean(RRpkWid(:,i))],'Color','k','LineWidth',2);
   scatter(zeros(numFlies,1)+20*i,40*squeeze(LRpkWid(:,i)),'filled','c'); 
   line([20*i-2 20*i+2],[40*mean(LRpkWid(:,i)) 40*mean(LRpkWid(:,i))],'Color','c','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 95]);
xlabel('vr bins (deg/sec)');
ylim([0 250]);

%% Align the peaks to one another as sorted by side and velocity - combining dark, gainOne, and gainTwo data
vRMin = 0;
vRMax = 360;
vRSpan = 30;
[LPB.dark.dataR, LPB.dark.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'dark');
[LPB.gainOne.dataR, LPB.gainOne.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'gainOne');
[LPB.gainTwo.dataR, LPB.gainTwo.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 0, 'gainTwo');
[RPB.dark.dataR, RPB.dark.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'dark');
[RPB.gainOne.dataR, RPB.gainOne.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'gainOne');
[RPB.gainTwo.dataR, RPB.gainTwo.dataG] = BumpAlign(allFlyData, vRMin, vRMax, vRSpan, 0, 1, 'gainTwo');

vRBinNum = round((vRMax-vRMin)/vRSpan);

PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID = 1:length(allFlyData)
    
    LPBDdataRAllStop = horzcat(horzcat(LPB.dark.dataR{flyID}.Stop,LPB.gainOne.dataR{flyID}.Stop),LPB.gainTwo.dataR{flyID}.Stop);
    LPBDdataGAllStop = horzcat(horzcat(LPB.dark.dataG{flyID}.Stop,LPB.gainOne.dataG{flyID}.Stop),LPB.gainTwo.dataG{flyID}.Stop);
    % Plot the dark data
    subplot(8,5,1+2*5*(flyID-1))
    hold on;
    plot(1:9,circshift(median(LPBDdataRAllStop,2),4),'r','LineWidth',2);
    plot(1:9,circshift(quantile(LPBDdataRAllStop,0.25,2),4),'r','LineWidth',1);
    plot(1:9,circshift(quantile(LPBDdataRAllStop,0.75,2),4),'r','LineWidth',1);
    plot(1:9,circshift(median(LPBDdataGAllStop,2),4),'g','LineWidth',2);
    plot(1:9,circshift(quantile(LPBDdataGAllStop,0.25,2),4),'g','LineWidth',1);
    plot(1:9,circshift(quantile(LPBDdataGAllStop,0.75,2),4),'g','LineWidth',1);
    ylim([1 2]);
    title(strcat('Left PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:4
        LPBDdataRAllCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CW{binID},LPB.gainOne.dataR{flyID}.CW{binID}),LPB.gainTwo.dataR{flyID}.CW{binID});
        LPBDdataGAllCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CW{binID},LPB.gainOne.dataG{flyID}.CW{binID}),LPB.gainTwo.dataG{flyID}.CW{binID});
        LPBDdataRAllCCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CCW{binID},LPB.gainOne.dataR{flyID}.CCW{binID}),LPB.gainTwo.dataR{flyID}.CCW{binID});
        LPBDdataGAllCCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CCW{binID},LPB.gainOne.dataG{flyID}.CCW{binID}),LPB.gainTwo.dataG{flyID}.CCW{binID});
        subplot(8,5,binID+1+2*5*(flyID-1));
        hold on;
        plot(1:9,circshift(median(LPBDdataRAllCW,2),4),'r','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataRAllCW,0.25,2),4),'r','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataRAllCW,0.75,2),4),'r','LineWidth',1);
        plot(1:9,circshift(median(LPBDdataGAllCW,2),4),'g','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataGAllCW,0.25,2),4),'g','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataGAllCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
        subplot(8,5,binID+6+2*5*(flyID-1));
        hold on;
        plot(1:9,circshift(median(LPBDdataRAllCCW,2),4),'r','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataRAllCCW,0.25,2),4),'r','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataRAllCCW,0.75,2),4),'r','LineWidth',1);
        plot(1:9,circshift(median(LPBDdataGAllCCW,2),4),'g','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataGAllCCW,0.25,2),4),'g','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataGAllCCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    
end
set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPosAll_Left'),'-dpdf');


PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID = 1:length(allFlyData)
    
    RPBDdataRAllStop = horzcat(horzcat(RPB.dark.dataR{flyID}.Stop,RPB.gainOne.dataR{flyID}.Stop),RPB.gainTwo.dataR{flyID}.Stop);
    RPBDdataGAllStop = horzcat(horzcat(RPB.dark.dataG{flyID}.Stop,RPB.gainOne.dataG{flyID}.Stop),RPB.gainTwo.dataG{flyID}.Stop);
    % Plot the dark data
    subplot(8,5,1+2*5*(flyID-1))
    hold on;
    plot(10:18,circshift(median(RPBDdataRAllStop,2),4),'r','LineWidth',2);
    plot(10:18,circshift(quantile(RPBDdataRAllStop,0.25,2),4),'r','LineWidth',1);
    plot(10:18,circshift(quantile(RPBDdataRAllStop,0.75,2),4),'r','LineWidth',1);
    plot(10:18,circshift(median(RPBDdataGAllStop,2),4),'g','LineWidth',2);
    plot(10:18,circshift(quantile(RPBDdataGAllStop,0.25,2),4),'g','LineWidth',1);
    plot(10:18,circshift(quantile(RPBDdataGAllStop,0.75,2),4),'g','LineWidth',1);
    ylim([1 2]);
    title(strcat('Left PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:4
        RPBDdataRAllCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CW{binID},RPB.gainOne.dataR{flyID}.CW{binID}),RPB.gainTwo.dataR{flyID}.CW{binID});
        RPBDdataGAllCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CW{binID},RPB.gainOne.dataG{flyID}.CW{binID}),RPB.gainTwo.dataG{flyID}.CW{binID});
        RPBDdataRAllCCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CCW{binID},RPB.gainOne.dataR{flyID}.CCW{binID}),RPB.gainTwo.dataR{flyID}.CCW{binID});
        RPBDdataGAllCCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CCW{binID},RPB.gainOne.dataG{flyID}.CCW{binID}),RPB.gainTwo.dataG{flyID}.CCW{binID});
        subplot(8,5,binID+1+2*5*(flyID-1));
        hold on;
        plot(10:18,circshift(median(RPBDdataRAllCW,2),4),'r','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataRAllCW,0.25,2),4),'r','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataRAllCW,0.75,2),4),'r','LineWidth',1);
        plot(10:18,circshift(median(RPBDdataGAllCW,2),4),'g','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataGAllCW,0.25,2),4),'g','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataGAllCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
        subplot(8,5,binID+6+2*5*(flyID-1));
        hold on;
        plot(10:18,circshift(median(RPBDdataRAllCCW,2),4),'r','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataRAllCCW,0.25,2),4),'r','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataRAllCCW,0.75,2),4),'r','LineWidth',1);
        plot(10:18,circshift(median(RPBDdataGAllCCW,2),4),'g','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataGAllCCW,0.25,2),4),'g','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataGAllCCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
end
set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPosAll_Right'),'-dpdf');

%% Align the peaks to one another as sorted by side and velocity - combining dark, gainOne, and gainTwo data and using only one peak to align
vRMin = 0;
vRMax = 360;
vRSpan = 30;
glomShift = 5;
[LPB.dark.dataR, RPB.dark.dataR, LPB.dark.dataG, RPB.dark.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'dark',glomShift);
[LPB.gainOne.dataR, RPB.gainOne.dataR, LPB.gainOne.dataG, RPB.gainOne.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainOne',glomShift);
[LPB.gainTwo.dataR, RPB.gainTwo.dataR, LPB.gainTwo.dataG, RPB.gainTwo.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainTwo',glomShift);

vRBinNum = round((vRMax-vRMin)/vRSpan);

PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID = 1:length(allFlyData)
    
    LPBDdataRAllStop = horzcat(horzcat(LPB.dark.dataR{flyID}.Stop,LPB.gainOne.dataR{flyID}.Stop),LPB.gainTwo.dataR{flyID}.Stop);
    LPBDdataGAllStop = horzcat(horzcat(LPB.dark.dataG{flyID}.Stop,LPB.gainOne.dataG{flyID}.Stop),LPB.gainTwo.dataG{flyID}.Stop);
    % Plot the dark data
    subplot(8,5,1+2*5*(flyID-1))
    hold on;
    plot(1:9,circshift(median(LPBDdataRAllStop,2),4),'r','LineWidth',2);
    plot(1:9,circshift(quantile(LPBDdataRAllStop,0.25,2),4),'r','LineWidth',1);
    plot(1:9,circshift(quantile(LPBDdataRAllStop,0.75,2),4),'r','LineWidth',1);
    plot(1:9,circshift(median(LPBDdataGAllStop,2),4),'g','LineWidth',2);
    plot(1:9,circshift(quantile(LPBDdataGAllStop,0.25,2),4),'g','LineWidth',1);
    plot(1:9,circshift(quantile(LPBDdataGAllStop,0.75,2),4),'g','LineWidth',1);
    ylim([1 2]);
    title(strcat('Left PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:4
        LPBDdataRAllCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CW{binID},LPB.gainOne.dataR{flyID}.CW{binID}),LPB.gainTwo.dataR{flyID}.CW{binID});
        LPBDdataGAllCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CW{binID},LPB.gainOne.dataG{flyID}.CW{binID}),LPB.gainTwo.dataG{flyID}.CW{binID});
        LPBDdataRAllCCW = horzcat(horzcat(LPB.dark.dataR{flyID}.CCW{binID},LPB.gainOne.dataR{flyID}.CCW{binID}),LPB.gainTwo.dataR{flyID}.CCW{binID});
        LPBDdataGAllCCW = horzcat(horzcat(LPB.dark.dataG{flyID}.CCW{binID},LPB.gainOne.dataG{flyID}.CCW{binID}),LPB.gainTwo.dataG{flyID}.CCW{binID});
        subplot(8,5,binID+1+2*5*(flyID-1));
        hold on;
        plot(1:9,circshift(median(LPBDdataRAllCW,2),4),'r','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataRAllCW,0.25,2),4),'r','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataRAllCW,0.75,2),4),'r','LineWidth',1);
        plot(1:9,circshift(median(LPBDdataGAllCW,2),4),'g','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataGAllCW,0.25,2),4),'g','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataGAllCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
        subplot(8,5,binID+6+2*5*(flyID-1));
        hold on;
        plot(1:9,circshift(median(LPBDdataRAllCCW,2),4),'r','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataRAllCCW,0.25,2),4),'r','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataRAllCCW,0.75,2),4),'r','LineWidth',1);
        plot(1:9,circshift(median(LPBDdataGAllCCW,2),4),'g','LineWidth',2);
        plot(1:9,circshift(quantile(LPBDdataGAllCCW,0.25,2),4),'g','LineWidth',1);
        plot(1:9,circshift(quantile(LPBDdataGAllCCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
    
end
set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPosAll_Left_OnePeak'),'-dpdf');


PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID = 1:length(allFlyData)
    
    RPBDdataRAllStop = horzcat(horzcat(RPB.dark.dataR{flyID}.Stop,RPB.gainOne.dataR{flyID}.Stop),RPB.gainTwo.dataR{flyID}.Stop);
    RPBDdataGAllStop = horzcat(horzcat(RPB.dark.dataG{flyID}.Stop,RPB.gainOne.dataG{flyID}.Stop),RPB.gainTwo.dataG{flyID}.Stop);
    % Plot the dark data
    subplot(8,5,1+2*5*(flyID-1))
    hold on;
    plot(10:18,circshift(median(RPBDdataRAllStop,2),4),'r','LineWidth',2);
    plot(10:18,circshift(quantile(RPBDdataRAllStop,0.25,2),4),'r','LineWidth',1);
    plot(10:18,circshift(quantile(RPBDdataRAllStop,0.75,2),4),'r','LineWidth',1);
    plot(10:18,circshift(median(RPBDdataGAllStop,2),4),'g','LineWidth',2);
    plot(10:18,circshift(quantile(RPBDdataGAllStop,0.25,2),4),'g','LineWidth',1);
    plot(10:18,circshift(quantile(RPBDdataGAllStop,0.75,2),4),'g','LineWidth',1);
    ylim([1 2]);
    title(strcat('Right PB, Fly ',num2str(flyID)));
    ylabel('Dark');
    xlabel('0');
    for binID=1:4
        RPBDdataRAllCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CW{binID},RPB.gainOne.dataR{flyID}.CW{binID}),RPB.gainTwo.dataR{flyID}.CW{binID});
        RPBDdataGAllCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CW{binID},RPB.gainOne.dataG{flyID}.CW{binID}),RPB.gainTwo.dataG{flyID}.CW{binID});
        RPBDdataRAllCCW = horzcat(horzcat(RPB.dark.dataR{flyID}.CCW{binID},RPB.gainOne.dataR{flyID}.CCW{binID}),RPB.gainTwo.dataR{flyID}.CCW{binID});
        RPBDdataGAllCCW = horzcat(horzcat(RPB.dark.dataG{flyID}.CCW{binID},RPB.gainOne.dataG{flyID}.CCW{binID}),RPB.gainTwo.dataG{flyID}.CCW{binID});
        subplot(8,5,binID+1+2*5*(flyID-1));
        hold on;
        plot(10:18,circshift(median(RPBDdataRAllCW,2),4),'r','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataRAllCW,0.25,2),4),'r','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataRAllCW,0.75,2),4),'r','LineWidth',1);
        plot(10:18,circshift(median(RPBDdataGAllCW,2),4),'g','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataGAllCW,0.25,2),4),'g','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataGAllCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
        subplot(8,5,binID+6+2*5*(flyID-1));
        hold on;
        plot(10:18,circshift(median(RPBDdataRAllCCW,2),4),'r','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataRAllCCW,0.25,2),4),'r','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataRAllCCW,0.75,2),4),'r','LineWidth',1);
        plot(10:18,circshift(median(RPBDdataGAllCCW,2),4),'g','LineWidth',2);
        plot(10:18,circshift(quantile(RPBDdataGAllCCW,0.25,2),4),'g','LineWidth',1);
        plot(10:18,circshift(quantile(RPBDdataGAllCCW,0.75,2),4),'g','LineWidth',1);
        ylim([1 2]);
    end
end
set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PBFig,strcat(allFlyData{flyID}.ID,'_vRotVsPBPosAll_Right_OnePeak'),'-dpdf');

%% Peak shift sanity check

flyID = 2;
trial = 'dark';
trialID = 3;
tSteps = [10:10:60];

RSig = allFlyData{flyID}.(trial){trialID}.RROIaveMax;
GSig = allFlyData{flyID}.(trial){trialID}.GROIaveMax;

figure;
for tStep = tSteps;
    subplot(2,length(tSteps),find(tSteps == tStep));
    hold on;
    plot(RSig(:,tStep+1),'r');
    plot(GSig(:,tStep+1),'g');
    ylim([1 2]);
    
    pkMax = max(RSig(1:8,tStep+1));
    pkShift = -find(RSig(1:9,tStep+1) == pkMax)+1;

    RSigNowLPB = circshift(RSig(1:9,tStep+1),pkShift);
    RSigNowRPB = circshift(RSig(10:18,tStep+1),pkShift);
    GSigNowLPB = circshift(GSig(1:9,tStep+1),pkShift);
    GSigNowRPB = circshift(GSig(10:18,tStep+1),pkShift);
    
    subplot(2,length(tSteps),length(tSteps)+find(tSteps == tStep));
    hold on;
    plot(vertcat(RSigNowLPB,RSigNowRPB),'r');
    plot(vertcat(GSigNowLPB,GSigNowRPB),'g');
    ylim([1 2]);
    
end

%% Find and plot the bump widths
glomShift = 5;
vRMin = 0;
vRMax = 720;
vRSpan = 30;
[LPB.dark.dataR, RPB.dark.dataR, LPB.dark.dataG, RPB.dark.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'dark',glomShift);
[LPB.gainOne.dataR, RPB.gainOne.dataR, LPB.gainOne.dataG, RPB.gainOne.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainOne',glomShift);
[LPB.gainTwo.dataR, RPB.gainTwo.dataR, LPB.gainTwo.dataG, RPB.gainTwo.dataG] = BumpAlignOnePeak(allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'gainTwo',glomShift);

vRBinNum = round((vRMax-vRMin)/vRSpan);

num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

LRpkWidStop = [];
LRpkWidStopSD = [];
RRpkWidStop = [];
RRpkWidStopSD = [];
LGpkWidStop = [];
LGpkWidStopSD = [];
RGpkWidStop = [];
RGpkWidStopSD = [];

LRpkWidCW = [];
LRpkWidCWSD = [];
RRpkWidCW = [];
RRpkWidCWSD = [];
LGpkWidCW = [];
LGpkWidCWSD = [];
RGpkWidCW = [];
RGpkWidCWSD = [];

LRpkWidCCW = [];
LRpkWidCCWSD = [];
RRpkWidCCW = [];
RRpkWidCCWSD = [];
LGpkWidCCW = [];
LGpkWidCCWSD = [];
RGpkWidCCW = [];
RGpkWidCCWSD = [];

SDFilt = 1;

plotRange = 270/vRSpan;
for flyID = 1:length(allFlyData)
    
    GallDat = [];
    RallDat = [];
    for trial = 1:length(allFlyData{flyID}.dark)
        GallDat = horzcat(GallDat,...
            reshape(allFlyData{flyID}.dark{trial}.GROIaveMax(greenSpan,:),...
            [1, length(greenSpan)*...
            size(allFlyData{flyID}.dark{trial}.GROIaveMax,2)]));
        RallDat = horzcat(RallDat,...
            reshape(allFlyData{flyID}.dark{trial}.RROIaveMax(redSpan,:),...
            [1,  length(redSpan)*...
            size(allFlyData{flyID}.dark{trial}.RROIaveMax,2)]));
    end
    GSDNow = std(GallDat);
    RSDNow = std(RallDat);
    
    LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
    RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
    LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
    RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;
    
    LRpkWids = zeros(size(LPBDataRAllStop,2),1);
    RRpkWids = zeros(size(RPBDataRAllStop,2),1);
    LGpkWids = zeros(size(LPBDataGAllStop,2),1);
    RGpkWids = zeros(size(RPBDataGAllStop,2),1);
    for StopIt = 1:size(LPBDataRAllStop,2)
       LRXCNow = LPBDataRAllStop(redSpan(1:8),StopIt);
       RRXCNow = RPBDataRAllStop(greenSpan(1:8),StopIt);
       LGXCNow = LPBDataGAllStop(greenSpan(1:8),StopIt);
       RGXCNow = RPBDataGAllStop(redSpan(1:8),StopIt);
       if max(LRXCNow) > mean(LRXCNow) + SDFilt*RSDNow
           overThresh = find((LRXCNow - min(LRXCNow))> 0.5*(max(LRXCNow)-min(LRXCNow)));
           LRpkWids(StopIt) = length(overThresh)-1;
       end
       if max(RRXCNow) > mean(RRXCNow) + SDFilt*RSDNow
           overThresh = find((RRXCNow - min(RRXCNow))> 0.5*(max(RRXCNow)-min(RRXCNow)));
           RRpkWids(StopIt) = length(overThresh)-1;
       end
       if max(LGXCNow) > mean(LGXCNow) + SDFilt*GSDNow
           overThresh = find((LGXCNow - min(LGXCNow))> 0.5*(max(LGXCNow)-min(LGXCNow)));
           LGpkWids(StopIt) = length(overThresh)-1;
       end
       if max(RGXCNow) > mean(RGXCNow) + SDFilt*GSDNow
           overThresh = find((RGXCNow - min(RGXCNow))> 0.5*(max(RGXCNow)-min(RGXCNow)));
           RGpkWids(StopIt) = length(overThresh)-1;
       end
    end
    LRpkWidStop(flyID) = mean(LRpkWids(find(LRpkWids>0)));
    LRpkWidStopSD(flyID) = std(LRpkWids(find(LRpkWids>0)));
    RRpkWidStop(flyID) = mean(RRpkWids(find(RRpkWids>0)));
    RRpkWidStopSD(flyID) = std(RRpkWids(find(RRpkWids>0)));
    LGpkWidStop(flyID) = mean(LGpkWids(find(LGpkWids>0)));
    LGpkWidStopSD(flyID) = std(LGpkWids(find(LGpkWids>0)));
    RGpkWidStop(flyID) = mean(RGpkWids(find(RGpkWids>0)));
    RGpkWidStopSD(flyID) = std(RGpkWids(find(RGpkWids>0)));
    
    % Find the dark data widths
        
    for binID=1:plotRange
        LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
        RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
        LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
        RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};
        
        if ~isempty(LPBDataRAllCW)
            LRpkWids = zeros(size(LPBDataRAllCW,2),1);
            RRpkWids = zeros(size(RPBDataRAllCW,2),1);
            LGpkWids = zeros(size(LPBDataGAllCW,2),1);
            RGpkWids = zeros(size(RPBDataGAllCW,2),1);
            for CWIt = 1:size(LPBDataRAllCW,2)
               LRXCNow = LPBDataRAllCW(redSpan(1:8),CWIt);
               RRXCNow = RPBDataRAllCW(greenSpan(1:8),CWIt);
               LGXCNow = LPBDataGAllCW(greenSpan(1:8),CWIt);
               RGXCNow = RPBDataGAllCW(redSpan(1:8),CWIt);
               if max(LRXCNow) > mean(LRXCNow) + SDFilt*RSDNow
                   overThresh = find((LRXCNow - min(LRXCNow))> 0.5*(max(LRXCNow)-min(LRXCNow)));
                   LRpkWids(CWIt) = length(overThresh)-1;
               end
               if max(RRXCNow) > mean(RRXCNow) + SDFilt*RSDNow
                   overThresh = find((RRXCNow - min(RRXCNow))> 0.5*(max(RRXCNow)-min(RRXCNow)));
                   RRpkWids(CWIt) = length(overThresh)-1;
               end
               if max(LGXCNow) > mean(LGXCNow) + SDFilt*GSDNow
                   overThresh = find((LGXCNow - min(LGXCNow))> 0.5*(max(LGXCNow)-min(LGXCNow)));
                   LGpkWids(CWIt) = length(overThresh)-1;
               end
               if max(RGXCNow) > mean(RGXCNow) + SDFilt*GSDNow
                   overThresh = find((RGXCNow - min(RGXCNow))> 0.5*(max(RGXCNow)-min(RGXCNow)));
                   RGpkWids(CWIt) = length(overThresh)-1;
               end
            end
            LRpkWidCW(flyID,binID) = mean(LRpkWids(find(LRpkWids>0)));
            LRpkWidCWSD(flyID,binID) = std(LRpkWids(find(LRpkWids>0)));
            RRpkWidCW(flyID,binID) = mean(RRpkWids(find(RRpkWids>0)));
            RRpkWidCWSD(flyID,binID) = std(RRpkWids(find(RRpkWids>0)));
            LGpkWidCW(flyID,binID) = mean(LGpkWids(find(LGpkWids>0)));
            LGpkWidCWSD(flyID,binID) = std(LGpkWids(find(LGpkWids>0)));
            RGpkWidCW(flyID,binID) = mean(RGpkWids(find(RGpkWids>0)));
            RGpkWidCWSD(flyID,binID) = std(RGpkWids(find(RGpkWids>0)));
        end
        
        LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
        RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
        LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
        RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};
        
        if ~isempty(LPBDataRAllCCW)
            LRpkWids = zeros(size(LPBDataRAllCCW,2),1);
            RRpkWids = zeros(size(RPBDataRAllCCW,2),1);
            LGpkWids = zeros(size(LPBDataGAllCCW,2),1);
            RGpkWids = zeros(size(RPBDataGAllCCW,2),1);
            for CCWIt = 1:size(LPBDataRAllCCW,2)
               LRXCNow = LPBDataRAllCCW(redSpan(1:8),CCWIt);
               RRXCNow = RPBDataRAllCCW(greenSpan(1:8),CCWIt);
               LGXCNow = LPBDataGAllCCW(greenSpan(1:8),CCWIt);
               RGXCNow = RPBDataGAllCCW(redSpan(1:8),CCWIt);
               if max(LRXCNow) > mean(LRXCNow) + SDFilt*RSDNow
                   overThresh = find((LRXCNow - min(LRXCNow))> 0.5*(max(LRXCNow)-min(LRXCNow)));
                   LRpkWids(CCWIt) = length(overThresh)-1;
               end
               if max(RRXCNow) > mean(RRXCNow) + SDFilt*RSDNow
                   overThresh = find((RRXCNow - min(RRXCNow))> 0.5*(max(RRXCNow)-min(RRXCNow)));
                   RRpkWids(CCWIt) = length(overThresh)-1;
               end
               if max(LGXCNow) > mean(LGXCNow) + SDFilt*GSDNow
                   overThresh = find((LGXCNow - min(LGXCNow))> 0.5*(max(LGXCNow)-min(LGXCNow)));
                   LGpkWids(CCWIt) = length(overThresh)-1;
               end
               if max(RGXCNow) > mean(RGXCNow) + SDFilt*GSDNow
                   overThresh = find((RGXCNow - min(RGXCNow))> 0.5*(max(RGXCNow)-min(RGXCNow)));
                   RGpkWids(CCWIt) = length(overThresh)-1;
               end
            end
            LRpkWidCCW(flyID,binID) = mean(LRpkWids(find(LRpkWids>0)));
            LRpkWidCCWSD(flyID,binID) = std(LRpkWids(find(LRpkWids>0)));
            RRpkWidCCW(flyID,binID) = mean(RRpkWids(find(RRpkWids>0)));
            RRpkWidCCWSD(flyID,binID) = std(RRpkWids(find(RRpkWids>0)));
            LGpkWidCCW(flyID,binID) = mean(LGpkWids(find(LGpkWids>0)));
            LGpkWidCCWSD(flyID,binID) = std(LGpkWids(find(LGpkWids>0)));
            RGpkWidCCW(flyID,binID) = mean(RGpkWids(find(RGpkWids>0)));
            RGpkWidCCWSD(flyID,binID) = std(RGpkWids(find(RGpkWids>0)));
        end
    end
end

save(strcat('D:\Imaging\2Color\Results\PBWidths_G',greenLine,'_R',redLine),...
    'LRpkWidStop','LRpkWidStopSD','RRpkWidStop','RRpkWidStopSD','LGpkWidStop','LGpkWidStopSD','RGpkWidStop','RGpkWidStopSD',...
    'LRpkWidCW','LRpkWidCWSD','RRpkWidCW','RRpkWidCWSD','LGpkWidCW','LGpkWidCWSD','RGpkWidCW','RGpkWidCWSD',...
    'LRpkWidCCW','LRpkWidCCWSD','RRpkWidCCW','RRpkWidCCWSD','LGpkWidCCW','LGpkWidCCWSD','RGpkWidCCW','RGpkWidCCWSD');


PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
hold on;
scatter(zeros(length(LRpkWidStop),1),LRpkWidStop,'k','filled');
for fly = 1:length(LRpkWidStop)
    line([0 0],[LRpkWidStop(fly)-LRpkWidStopSD(fly), LRpkWidStop(fly)+LRpkWidStopSD(fly)],'Color','k');
end
for binID = 1:plotRange
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCCW(:,binID),'m','filled');
    for fly = 1:length(LRpkWidStop)
        line([binID binID],[LRpkWidCW(fly,binID)-LRpkWidCWSD(fly,binID), LRpkWidCW(fly,binID)+LRpkWidCWSD(fly,binID)],'Color','c');
        line([binID binID],[LRpkWidCCW(fly,binID)-LRpkWidCCWSD(fly,binID), LRpkWidCCW(fly,binID)+LRpkWidCCWSD(fly,binID)],'Color','m');
    end
end
ylim([0 3.5]);

subplot(2,2,2);
hold on;
scatter(zeros(length(RRpkWidStop),1),RRpkWidStop,'k','filled');
for fly = 1:length(RRpkWidStop)
    line([0 0],[RRpkWidStop(fly)-RRpkWidStopSD(fly), RRpkWidStop(fly)+RRpkWidStopSD(fly)],'Color','k');
end
for binID = 1:plotRange
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCCW(:,binID),'m','filled');
    for fly = 1:length(RRpkWidStop)
        line([binID binID],[RRpkWidCW(fly,binID)-RRpkWidCWSD(fly,binID), RRpkWidCW(fly,binID)+RRpkWidCWSD(fly,binID)],'Color','c');
        line([binID binID],[RRpkWidCCW(fly,binID)-RRpkWidCCWSD(fly,binID), RRpkWidCCW(fly,binID)+RRpkWidCCWSD(fly,binID)],'Color','m');
    end
end
ylim([0 3.5]);

subplot(2,2,3);
hold on;
scatter(zeros(length(LGpkWidStop),1),LGpkWidStop,'k','filled');
for fly = 1:length(LGpkWidStop)
    line([0 0],[LGpkWidStop(fly)-LGpkWidStopSD(fly), LGpkWidStop(fly)+LGpkWidStopSD(fly)],'Color','k');
end
for binID = 1:plotRange
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCCW(:,binID),'m','filled');
    for fly = 1:length(LRpkWidStop)
        line([binID binID],[LGpkWidCW(fly,binID)-LGpkWidCWSD(fly,binID), LGpkWidCW(fly,binID)+LGpkWidCWSD(fly,binID)],'Color','c');
        line([binID binID],[LGpkWidCCW(fly,binID)-LGpkWidCCWSD(fly,binID), LGpkWidCCW(fly,binID)+LGpkWidCCWSD(fly,binID)],'Color','m');
    end
end
ylim([0 3.5]);

subplot(2,2,4);
hold on;
scatter(zeros(length(RGpkWidStop),1),RGpkWidStop,'k','filled');
for fly = 1:length(RGpkWidStop)
    line([0 0],[RGpkWidStop(fly)-RGpkWidStopSD(fly), RGpkWidStop(fly)+RGpkWidStopSD(fly)],'Color','k');
end
for binID = 1:plotRange
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCCW(:,binID),'m','filled');
    for fly = 1:length(RRpkWidStop)
        line([binID binID],[RGpkWidCW(fly,binID)-RGpkWidCWSD(fly,binID), RGpkWidCW(fly,binID)+RGpkWidCWSD(fly,binID)],'Color','c');
        line([binID binID],[RGpkWidCCW(fly,binID)-RGpkWidCCWSD(fly,binID), RGpkWidCCW(fly,binID)+RGpkWidCCWSD(fly,binID)],'Color','m');
    end
end
ylim([0 3.5]);

%% Plot bump widths across color directions

PBFig = figure;

vRSpan = 30;
plotRange = 150/vRSpan;

load('D:\Imaging\2Color\Results\PBWidths_G60D05_R37F06.mat');

L60D05Stop = LGpkWidStop;
R60D05Stop = RGpkWidStop;
L37F06Stop = LRpkWidStop;
R37F06Stop = RRpkWidStop;

L60D05CW = LGpkWidCW;
R60D05CW = RGpkWidCW;
L37F06CW = LRpkWidCW;
R37F06CW = RRpkWidCW;

L60D05CCW = LGpkWidCCW;
R60D05CCW = RGpkWidCCW;
L37F06CCW = LRpkWidCCW;
R37F06CCW = RRpkWidCCW;

subplot(2,2,1);
hold on;
scatter(zeros(length(LRpkWidStop),1),LRpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(LRpkWidStop)
   line([0 1],[LRpkWidStop(fly) LRpkWidCW(fly,1)],'LineWidth',1,'Color','c');
   line([0 1],[LRpkWidStop(fly) LRpkWidCCW(fly,1)],'LineWidth',1,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[LRpkWidCW(fly,binID) LRpkWidCW(fly,binID+1)],'LineWidth',1,'Color','c');
       line([binID binID+1],[LRpkWidCCW(fly,binID) LRpkWidCCW(fly,binID+1)],'LineWidth',1,'Color','m');
   end
end

subplot(2,2,2);
hold on;
scatter(zeros(length(RRpkWidStop),1),RRpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(RRpkWidStop)
   line([0 1],[RRpkWidStop(fly) RRpkWidCW(fly,1)],'LineWidth',1,'Color','c');
   line([0 1],[RRpkWidStop(fly) RRpkWidCCW(fly,1)],'LineWidth',1,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[RRpkWidCW(fly,binID) RRpkWidCW(fly,binID+1)],'LineWidth',1,'Color','c');
       line([binID binID+1],[RRpkWidCCW(fly,binID) RRpkWidCCW(fly,binID+1)],'LineWidth',1,'Color','m');
   end
end

subplot(2,2,3);
hold on;
scatter(zeros(length(LGpkWidStop),1),LGpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(LGpkWidStop)
   line([0 1],[LGpkWidStop(fly) LGpkWidCW(fly,1)],'LineWidth',1,'Color','c');
   line([0 1],[LGpkWidStop(fly) LGpkWidCCW(fly,1)],'LineWidth',1,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[LGpkWidCW(fly,binID) LGpkWidCW(fly,binID+1)],'LineWidth',1,'Color','c');
       line([binID binID+1],[LGpkWidCCW(fly,binID) LGpkWidCCW(fly,binID+1)],'LineWidth',1,'Color','m');
   end
end

subplot(2,2,4);
hold on;
scatter(zeros(length(RGpkWidStop),1),RGpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(LGpkWidStop)
   line([0 1],[RGpkWidStop(fly) RGpkWidCW(fly,1)],'LineWidth',1,'Color','c');
   line([0 1],[RGpkWidStop(fly) RGpkWidCCW(fly,1)],'LineWidth',1,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[RGpkWidCW(fly,binID) RGpkWidCW(fly,binID+1)],'LineWidth',1,'Color','c');
       line([binID binID+1],[RGpkWidCCW(fly,binID) RGpkWidCCW(fly,binID+1)],'LineWidth',1,'Color','m');
   end
end


load('D:\Imaging\2Color\Results\PBWidths_G37F06_R60D05.mat');

L60D05Stop = horzcat(L60D05Stop,LRpkWidStop);
R60D05Stop = horzcat(R60D05Stop,RRpkWidStop);
L37F06Stop = horzcat(L37F06Stop,LGpkWidStop);
R37F06Stop = horzcat(R37F06Stop,RGpkWidStop);

L60D05CW = vertcat(L60D05CW,LRpkWidCW);
R60D05CW = vertcat(R60D05CW,RRpkWidCW);
L37F06CW = vertcat(L37F06CW,LGpkWidCW);
R37F06CW = vertcat(R37F06CW,RGpkWidCW);

L60D05CCW = vertcat(L60D05CCW,LRpkWidCCW);
R60D05CCW = vertcat(R60D05CCW,RRpkWidCCW);
L37F06CCW = vertcat(L37F06CCW,LGpkWidCCW);
R37F06CCW = vertcat(R37F06CCW,RGpkWidCCW);

subplot(2,2,1);
hold on;
scatter(zeros(length(LGpkWidStop),1),LGpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LGpkWidStop),1)+binID,LGpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(LGpkWidStop)
   line([0 1],[LGpkWidStop(fly) LGpkWidCW(fly,1)],'LineWidth',0.5,'Color','c');
   line([0 1],[LGpkWidStop(fly) LGpkWidCCW(fly,1)],'LineWidth',0.5,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[LGpkWidCW(fly,binID) LGpkWidCW(fly,binID+1)],'LineWidth',0.5,'Color','c');
       line([binID binID+1],[LGpkWidCCW(fly,binID) LGpkWidCCW(fly,binID+1)],'LineWidth',0.5,'Color','m');
   end
end

subplot(2,2,2);
hold on;
scatter(zeros(length(RGpkWidStop),1),RGpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RGpkWidStop),1)+binID,RGpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(RGpkWidStop)
   line([0 1],[RGpkWidStop(fly) RGpkWidCW(fly,1)],'LineWidth',0.5,'Color','c');
   line([0 1],[RGpkWidStop(fly) RGpkWidCCW(fly,1)],'LineWidth',0.5,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[RGpkWidCW(fly,binID) RGpkWidCW(fly,binID+1)],'LineWidth',0.5,'Color','c');
       line([binID binID+1],[RGpkWidCCW(fly,binID) RGpkWidCCW(fly,binID+1)],'LineWidth',0.5,'Color','m');
   end
end

subplot(2,2,3);
hold on;
scatter(zeros(length(LRpkWidStop),1),LRpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(LRpkWidStop),1)+binID,LRpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(LRpkWidStop)
   line([0 1],[LRpkWidStop(fly) LRpkWidCW(fly,1)],'LineWidth',0.5,'Color','c');
   line([0 1],[LRpkWidStop(fly) LRpkWidCCW(fly,1)],'LineWidth',0.5,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[LRpkWidCW(fly,binID) LRpkWidCW(fly,binID+1)],'LineWidth',0.5,'Color','c');
       line([binID binID+1],[LRpkWidCCW(fly,binID) LRpkWidCCW(fly,binID+1)],'LineWidth',0.5,'Color','m');
   end
end

subplot(2,2,4);
hold on;
scatter(zeros(length(RRpkWidStop),1),RRpkWidStop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCW(:,binID),'c','filled');
    scatter(zeros(length(RRpkWidStop),1)+binID,RRpkWidCCW(:,binID),'m','filled');
end
for fly = 1:length(RRpkWidStop)
   line([0 1],[RRpkWidStop(fly) RRpkWidCW(fly,1)],'LineWidth',0.5,'Color','c');
   line([0 1],[RRpkWidStop(fly) RRpkWidCCW(fly,1)],'LineWidth',0.5,'Color','m');
   for binID = 1:plotRange-1
       line([binID binID+1],[RRpkWidCW(fly,binID) RRpkWidCW(fly,binID+1)],'LineWidth',0.5,'Color','c');
       line([binID binID+1],[RRpkWidCCW(fly,binID) RRpkWidCCW(fly,binID+1)],'LineWidth',0.5,'Color','m');
   end
end

subplot(2,2,1);
line([-0.5 0.5], [mean(L37F06Stop(~isnan(L37F06Stop))) mean(L37F06Stop(~isnan(L37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    line([binID-0.5 binID+0.5], [mean(L37F06CW(:,binID)) mean(L37F06CW(:,binID))],'LineWidth',2,'Color','c');
    line([binID-0.5 binID+0.5], [mean(L37F06CCW(:,binID)) mean(L37F06CCW(:,binID))],'LineWidth',2,'Color','m');
end
xlim([-0.5 plotRange+0.5]);
ylim([0 3.5]);
ylabel('half width(# of glomeruli)');

subplot(2,2,2);
line([-0.5 0.5], [mean(R37F06Stop(~isnan(R37F06Stop))) mean(R37F06Stop(~isnan(R37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    line([binID-0.5 binID+0.5], [mean(R37F06CW(:,binID)) mean(R37F06CW(:,binID))],'LineWidth',2,'Color','c');
    line([binID-0.5 binID+0.5], [mean(R37F06CCW(:,binID)) mean(R37F06CCW(:,binID))],'LineWidth',2,'Color','m');
end
xlim([-0.5 plotRange+0.5]);
ylim([0 3.5]);

subplot(2,2,3);
line([-0.5 0.5], [mean(L60D05Stop(~isnan(L60D05Stop))) mean(L60D05Stop(~isnan(L60D05Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    line([binID-0.5 binID+0.5], [mean(L60D05CW(:,binID)) mean(L60D05CW(:,binID))],'LineWidth',2,'Color','c');
    line([binID-0.5 binID+0.5], [mean(L60D05CCW(:,binID)) mean(L60D05CCW(:,binID))],'LineWidth',2,'Color','m');
end
xlim([-0.5 plotRange+0.5]);
ylim([0 3.5]);
ylabel('half width(# of glomeruli)');

subplot(2,2,4);
line([-0.5 0.5], [mean(R60D05Stop(~isnan(R60D05Stop))) mean(R60D05Stop(~isnan(R60D05Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    line([binID-0.5 binID+0.5], [mean(R60D05CW(:,binID)) mean(R60D05CW(:,binID))],'LineWidth',2,'Color','c');
    line([binID-0.5 binID+0.5], [mean(R60D05CCW(:,binID)) mean(R60D05CCW(:,binID))],'LineWidth',2,'Color','m');
end
xlim([-0.5 plotRange+0.5]);
ylim([0 3.5]);

%% Plot the bump widths across color directions - consolidated L and R, CW and CCW
PBFig = figure;

vRSpan = 30;
plotRange = 150/vRSpan;

load('D:\Imaging\2Color\Results\PBWidths_G60D05_R37F06.mat');

allG60D05Stop = (LGpkWidStop+RGpkWidStop)/2;
allR37F06Stop = (LRpkWidStop+RRpkWidStop)/2;

allG60D05Mov = (LGpkWidCW+RGpkWidCW+LGpkWidCCW+RGpkWidCCW)/4;
allR37F06Mov = (LRpkWidCW+RRpkWidCW+LRpkWidCCW+RRpkWidCCW)/4;

subplot(2,2,1);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(allG60D05Stop),1)+binID,allG60D05Mov(:,binID),'g','filled');
end
for fly = 1:size(allG60D05Mov,1)
    plot(1:plotRange,allG60D05Mov(fly,1:plotRange),'color','g')
end


subplot(2,2,2);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(allR37F06Stop),1)+binID,allR37F06Mov(:,binID),'r','filled');
end
for fly = 1:size(allR37F06Mov,1)
    plot(1:plotRange,allR37F06Mov(fly,1:plotRange),'color','r')
end

load('D:\Imaging\2Color\Results\PBWidths_G37F06_R60D05.mat');

allR60D05Stop = (LRpkWidStop+RRpkWidStop)/2;
allG37F06Stop = (LGpkWidStop+RGpkWidStop)/2;

allR60D05Mov = (LRpkWidCW+RRpkWidCW+LRpkWidCCW+RRpkWidCCW)/4;
allG37F06Mov = (LGpkWidCW+RGpkWidCW+LGpkWidCCW+RGpkWidCCW)/4;

subplot(2,2,1);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(allR60D05Stop),1)+binID,allR60D05Mov(:,binID),'r','filled');
end
for fly = 1:size(allR60D05Mov,1)
    plot(1:plotRange,allR60D05Mov(fly,1:plotRange),'color','r')
end

subplot(2,2,2);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    scatter(zeros(length(allG37F06Stop),1)+binID,allG37F06Mov(:,binID),'g','filled');
end
for fly = 1:size(allG37F06Mov,1)
    plot(1:plotRange,allG37F06Mov(fly,1:plotRange),'color','g')
end

subplot(2,2,1);
% line([-0.5 0.5], [mean(L37F06Stop(~isnan(L37F06Stop))) mean(L37F06Stop(~isnan(L37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    all4Bin = vertcat(allR60D05Mov(:,binID),allG60D05Mov(:,binID));
    all4Bin = all4Bin(~isnan(all4Bin));
    line([binID-0.25 binID+0.25], [mean(all4Bin(find(all4Bin))) mean(all4Bin(find(all4Bin)))],....
        'LineWidth',2,'Color','k');
end
xlim([0.5 plotRange+0.5]);
ylim([0 3.5]);
ylabel('half width(# of glomeruli)');

subplot(2,2,2);
% line([-0.5 0.5], [mean(R37F06Stop(~isnan(R37F06Stop))) mean(R37F06Stop(~isnan(R37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    all4Bin = vertcat(allR37F06Mov(:,binID),allG37F06Mov(:,binID));
    all4Bin = all4Bin(~isnan(all4Bin));
    line([binID-0.25 binID+0.25], [mean(all4Bin(find(all4Bin))) mean(all4Bin(find(all4Bin)))],....
        'LineWidth',2,'Color','k');
end
xlim([0.5 plotRange+0.5]);
ylim([0 3.5]);

%% As above, but for VT8135
PBFig = figure;

vRSpan = 30;
plotRange = 150/vRSpan;

load('D:\Imaging\2Color\Results\PBWidths_G60D05_RVT8135.mat');

allG60D05Stop = (LGpkWidStop+RGpkWidStop)/2;
allR37F06Stop = (LRpkWidStop+RRpkWidStop)/2;

allG60D05Mov = (LGpkWidCW+RGpkWidCW+LGpkWidCCW+RGpkWidCCW)/4;
allR37F06Mov = (LRpkWidCW+RRpkWidCW+LRpkWidCCW+RRpkWidCCW)/4;

subplot(2,2,1);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    GNow = allG60D05Mov(:,binID);
    GNow = GNow(find(GNow));
    scatter(zeros(length(GNow),1)+binID,GNow,'g','filled');
end
for fly = 1:size(allG60D05Mov,1)
    plot(1:plotRange,allG60D05Mov(fly,1:plotRange),'color','g')
end


subplot(2,2,2);
hold on;
% scatter(zeros(length(all60D05Stop),1),all60D05Stop,'k','filled');
for binID = 1:plotRange
    RNow = allR37F06Mov(:,binID);
    RNow = RNow(find(RNow));
    scatter(zeros(length(RNow),1)+binID,RNow,'r','filled');
end
for fly = 1:size(allR37F06Mov,1)
    plot(1:plotRange,allR37F06Mov(fly,1:plotRange),'color','r')
end


subplot(2,2,1);
% line([-0.5 0.5], [mean(L37F06Stop(~isnan(L37F06Stop))) mean(L37F06Stop(~isnan(L37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    all4Bin = allG60D05Mov(:,binID);
    all4Bin = all4Bin(~isnan(all4Bin));
    all4Bin = all4Bin(find(all4Bin));
    line([binID-0.25 binID+0.25], [mean(all4Bin(find(all4Bin))) mean(all4Bin(find(all4Bin)))],....
        'LineWidth',2,'Color','k');
end
xlim([0.5 plotRange+0.5]);
ylim([0 3.5]);
ylabel('half width(# of glomeruli)');

subplot(2,2,2);
% line([-0.5 0.5], [mean(R37F06Stop(~isnan(R37F06Stop))) mean(R37F06Stop(~isnan(R37F06Stop)))],'LineWidth',2,'Color','k');
for binID = 1:plotRange
    all4Bin = allR37F06Mov(:,binID);
    all4Bin = all4Bin(~isnan(all4Bin));
    all4Bin = all4Bin(find(all4Bin));
    line([binID-0.25 binID+0.25], [mean(all4Bin(find(all4Bin))) mean(all4Bin(find(all4Bin)))],....
        'LineWidth',2,'Color','k');
end
xlim([0.5 plotRange+0.5]);
ylim([0 3.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to analyze all of the 2Color EB data (already processed with ProcessAllData2Color) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear out old data
clear;
clc;
    
%% Load all of the data
% Specify directories
dirRoot = 'D:\Imaging\2Color';
% dirRoot = 'D:\Imaging\2Color\LightBall';
% 
% imDir{1} = '20170112';
% imDir{2} = '20170114';
% redLine = '60D05';
% greenLine = '37F06';

% imDir{1} = '20160504';
% redLine = '37F06';
% greenLine = '60D05';
% 
imDir{1} = '20160423-Flip';
imDir{2} = '20160425-Flip';
imDir{3} = '20160426-Flip';
imDir{4} = '20160429';
imDir{5} = '20160505';
redLine = '60D05';
greenLine = '37F06';

% imDir{1} = '20160817';
% redLine = '60D05';
% greenLine = 'VT20739';

% imDir{1} = '20160618\EB';
% redLine = 'VT20739';
% greenLine = '60D05';

% imDir{1} = '20160815';
% redLine = 'VT45661';
% greenLine = '60D05';

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

% Get the position data, ROI data, and trial ID for all of the flies and all of the trials
% for dirs = 1:length(imDir)
%     fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
%     for fileID = 4:length(fileNames)
%         fileName = fileNames(fileID).name;
%         if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs')
%             flyName = strcat(imDir{dirs},'-',fileName(1:4));
%             load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
%             for flyID = 1:length(allFlyData)
%                 if strcmp(flyName,allFlyData{flyID}.ID)
%                     fileParts = strsplit(fileName,'_');
%                     if strcmpi(fileParts(5),'1x')
%                         allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.positionDat = positionDat;
%                         allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMax = RROIaveMax;
%                         allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMax = GROIaveMax;
%                         allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMean;
%                         allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMean;
%                     elseif strcmpi(fileParts(5),'2x')
%                         allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.positionDat = positionDat;
%                         allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMax = RROIaveMax;
%                         allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMax = GROIaveMax;
%                         allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMean = RROIaveMean;
%                         allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMean = GROIaveMean;
%                     elseif strcmpi(fileParts(5),'Dark')
%                         allFlyData{flyID}.dark{str2num(fileName(end-4))}.positionDat = positionDat;
%                         allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMax = RROIaveMax;
%                         allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMax = GROIaveMax;
%                         allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMean = RROIaveMean;
%                         allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMean = GROIaveMean;  
%                     end
%                 end
%             end
%         end 
%     end
% end

% Get the position data, ROI data, and trial ID for all of the flies and all of the trials
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
%             RROIaveMeanFilt = sgolayfilt(RROIaveMean,sgolayOrder,sgolayWindow,[],2);
%             GROIaveMeanFilt = sgolayfilt(GROIaveMean,sgolayOrder,sgolayWindow,[],2);
            
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
%                             allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
%                             allFlyData{flyID}.gainOne{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'2x')
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
%                             allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
%                             allFlyData{flyID}.gainTwo{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);
                        elseif strcmpi(fileParts(5),'Dark')
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.positionDatMatch = positionDatMatch;
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                            allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
%                             allFlyData{flyID}.dark{str2num(fileName(end-4))}.RROIaveMean = RROIaveMeanFilt(:,minFG:maxFG);
%                             allFlyData{flyID}.dark{str2num(fileName(end-4))}.GROIaveMean = GROIaveMeanFilt(:,minFG:maxFG);  
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

            % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(RROIaveMaxFilt);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(GROIaveMaxFilt);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(RROIaveMaxFilt) 3]);
            overlayIm(:,:,1) = (RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
            overlayIm(:,:,2) = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            subplot(4,1,4);
            plot(OffsetRotMatch(:,1),-OffsetRotMatch(:,2),'LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([0 OffsetRotMatch(end,1)]);
            ylim([-pi pi]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_Dark_Summary_',num2str(darkID)),'-dpdf');
            print(RGcomp,strcat('Dark_Summary_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
    
    % Plot parameters for the 1x gain data
    for gainOneID = 1:length(allFlyData{flyID}.gainOne)
        if ~isempty(allFlyData{flyID}.gainOne{gainOneID});
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax));
            tSpan = allFlyData{flyID}.gainOne{gainOneID}.positionDat.t - ...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(1) + ...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab)/num_planes);

            % Match the position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainOne{gainOneID}.positionDat.t >=...
                    (allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetLat(tMatch(1));
            end

            % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(RROIaveMaxFilt);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(GROIaveMaxFilt);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(RROIaveMaxFilt) 3]);
            overlayIm(:,:,1) = (RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
            overlayIm(:,:,2) = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            subplot(4,1,4);
            plot(OffsetRotMatch(:,1),-OffsetRotMatch(:,2),'LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([0 OffsetRotMatch(end,1)]);
            ylim([-pi pi]);

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_1x_Summary_',num2str(gainOneID)),'-dpdf');
            print(RGcomp,strcat('1x_Summary_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
    
    % Plot parameters for the 2x gain data
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo)
        if ~isempty(allFlyData{flyID}.gainTwo{gainTwoID});
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax));
            tSpan = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t - ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) + ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/num_planes);

            % Match the stripe position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t >=...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetLat(tMatch(1));
            end
            
            stripeMatch = mod(2*UnWrap(OffsetRotMatch(:,2),2,0),2*pi)-pi;
            
%             [posRot, posFor, posLat] = PositionConverter(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat);
%             
%             % Match the fly position data to the framegrab times
%             posRotMatch = zeros(maxFG-minFG+1,1);
%             posForMatch = zeros(maxFG-minFG+1,1);
%             posLatMatch = zeros(maxFG-minFG+1,1);
%             for interp = minFG:maxFG
%                 tMatch = find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t >=...
%                     (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) +...
%                     (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
%                 posRotMatch(interp-minFG+1,2) = pi/180*posRot(tMatch(1));
%                 posForMatch(interp-minFG+1) = posFor(tMatch(1));
%                 posLatMatch(interp-minFG+1) = posLat(tMatch(1));
%             end

            % Plot the relative values
            RGcomp = figure('units','normalized','outerposition',[0 0 1 1]);

            subplot(4,1,1);
            imagesc(RROIaveMaxFilt);
            axis off;
            title(redLine);

            subplot(4,1,2);
            imagesc(GROIaveMaxFilt);
            axis off;
            title(greenLine);

            subplot(4,1,3);
            overlayIm = zeros([size(RROIaveMaxFilt) 3]);
            overlayIm(:,:,1) = (RROIaveMaxFilt-min(min(RROIaveMaxFilt)))./(max(max(RROIaveMaxFilt))-min(min(RROIaveMaxFilt)));
            overlayIm(:,:,2) = (GROIaveMaxFilt-min(min(GROIaveMaxFilt)))./(max(max(GROIaveMaxFilt))-min(min(GROIaveMaxFilt)));
            image(overlayIm);
            axis off;

            colormap(brewermap(64, 'Blues'));

            subplot(4,1,4);
            hold on;
            plot(OffsetRotMatch(:,1),-OffsetRotMatch(:,2),'LineWidth',2,'LineStyle','--');
            plot(OffsetRotMatch(:,1),-stripeMatch,'LineWidth',2);
            axis tight;
            ylabel('Global rotation (rad)','FontSize',16);
            xlabel('Time (sec)','FontSize',16);
            set(gca,'FontSize',16);
            xlim([0 OffsetRotMatch(end,1)]);
            ylim([-pi pi]);
            legend('Stripe (x2)','Fly');

            set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_2x_Summary_',num2str(gainTwoID)),'-dpdf');
            print(RGcomp,strcat('2x_Summary_',num2str(darkID)),'-dpdf');
            
            delete(RGcomp);
            
            clear RROIaveMeanFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

%% Plot the flies' trajectories for the dark trials
flyTraj = figure('units','normalized','outerposition',[0 0 1 1]);

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    vf = [];
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID})
            subplot(numFlies,7,darkID+(flyID-1)*7);
            num_planes = round(length(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.dark{darkID}.RROIaveMax));
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
            plot(-OffsetLatMatch,OffsetForMatch,'Color',[1-flyID/numFlies flyID/numFlies, flyID/numFlies]);
            xlim([-5 5]);
            ylim([-5 5]);
            set(gca,'FontSize',12);
            xlabel('cm');
            
            trtsity = 0;
            ptsNCurve = round(4./mean(diff(OffsetRotMatch(:,1))));
            TransMatch = sqrt(diff(OffsetLatMatch).^2+diff(OffsetForMatch).^2);
            for i=1:length(OffsetLatMatch)-ptsNCurve
                L = sqrt((OffsetLatMatch(i)-OffsetLatMatch(i+ptsNCurve))^2 ...
                    +(OffsetForMatch(i)-OffsetForMatch(i+ptsNCurve))^2);
                if L > 0
                    trtsity = trtsity + abs(sum(TransMatch(i:i+ptsNCurve-1)))./L;
                end
            end
            trtsity = trtsity./(length(OffsetLatMatch)-ptsNCurve);
            text(2, 4, num2str(trtsity));

            if darkID == 1
                 ylabel(allFlyData{flyID}.ID);
            else
                ylabel('cm');
            end
        end
    end
end

set(flyTraj,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(flyTraj,strcat('Results\','FlyTraj_Red',redLine,'_Green',greenLine),'-dpdf');
            
%% Compare the bump shifts for different velocity thresholds for the dark trials
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
velThreshLC = [0 20 40 60];
velThreshUC = [20 40 60 80];
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
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(tSmooth./mean(diff(OffsetRotMatch(:,1)))));
            
            % Find the PVA for considering only periods where the signals
            % are strong
            % Calculate the PVA for 60D05
            num_ROIs = 16;
            angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
            angsraw = angsraw';
            clear meanAngRawG meanAngRawR;
            clear meanIntRawG meanIntRawR;
            for ts = 1:length(GROIaveMaxFilt)
                meanAngRawG(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
                meanAngRawR(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
                meanIntRawG(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
                meanIntRawR(ts) = circ_r(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
            end
                    
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins
                RturnsG = [];
                RturnsR = [];
                LturnsG = [];
                LturnsR = [];

                velThreshUCNow = velThreshUC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                velThreshLCNow = velThreshLC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));

                for i=minFG+1:maxFG
                    if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < velThreshUCNow
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        RturnsG = horzcat(RturnsG, GROIaveMaxFilt(:,i));
                        RturnsR = horzcat(RturnsR, RROIaveMaxFilt(:,i));
                    elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > -velThreshUCNow
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        LturnsG = horzcat(LturnsG, GROIaveMaxFilt(:,i));
                        LturnsR = horzcat(LturnsR, RROIaveMaxFilt(:,i));
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsG = RturnsG;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsR = RturnsR;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsG = LturnsG;
                XCVals{flyID}.vel{velCut}.dark{darkID}.LturnsR = LturnsR;
            end
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       RturnsGallRalign = [];
       RturnsGallGalign = [];
       RturnsRallRalign = [];
       RturnsRallGalign = [];
       LturnsGallRalign = [];
       LturnsGallGalign = [];
       LturnsRallRalign = [];
       LturnsRallGalign = [];
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
            RturnsGRalign = [];
            RturnsGGalign = [];
            RturnsRRalign = [];
            RturnsRGalign = [];
            LturnsGRalign = [];
            LturnsGGalign = [];
            LturnsRRalign = [];
            LturnsRGalign = [];
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.dark{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsG(:,i)));
                    RturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsG(:,i),num_ROIs/2-pkPosRalign);
                    RturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsG(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsR(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsR(:,i),num_ROIs/2-pkPosGalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsG(:,i)));
                    LturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsG(:,i),num_ROIs/2-pkPosRalign);
                    LturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsG(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsR(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsR(:,i),num_ROIs/2-pkPosGalign);
                end
            end
            
            RturnsGallRalign = horzcat(RturnsGallRalign,RturnsGRalign);
            RturnsGallGalign = horzcat(RturnsGallGalign,RturnsGGalign);
            RturnsRallRalign = horzcat(RturnsRallRalign,RturnsRRalign);
            RturnsRallGalign = horzcat(RturnsRallGalign,RturnsRGalign);
            LturnsGallRalign = horzcat(LturnsGallRalign,LturnsGRalign);
            LturnsGallGalign = horzcat(LturnsGallGalign,LturnsGGalign);
            LturnsRallRalign = horzcat(LturnsRallRalign,LturnsRRalign);
            LturnsRallGalign = horzcat(LturnsRallGalign,LturnsRGalign);

       end
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalign = RturnsGallRalign;
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalign = RturnsGallGalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalign = RturnsRallRalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalign = RturnsRallGalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalign = LturnsGallRalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalign = LturnsGallGalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalign = LturnsRallRalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalign = LturnsRallGalign;
    end 
end

bumpCompG = figure('units','normalized','outerposition',[0 0 0.5 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalign;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetG,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetG,0.75,2),'b','LineWidth',1);
           plot(median(RnetG,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetG,2)),'Color','b');
           plot(median(RnetR,2),'--r','LineWidth',1);
           
           plot(quantile(LnetG,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetG,0.75,2),'m','LineWidth',1);
           plot(median(LnetG,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetG,2)),'Color','m');
           plot(median(LnetR,2),'--r','LineWidth',1);
           
           xlim([1 16]);
           ylim([-0.5 2]);
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

bumpCompR = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalign;
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
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           plot(median(LnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([-0.5 2]);
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
% print(bumpCompR,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine),'-dpdf');
print(bumpCompR,strcat('BumpCompGAlign_Red',redLine,'_Green',greenLine),'-dpdf');
set(bumpCompG,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompG,strcat('Results\BumpCompRAlign_Red',redLine,'_Green',greenLine),'-dpdf');
print(bumpCompG,strcat('BumpCompRAlign_Red',redLine,'_Green',greenLine),'-dpdf');

%% Compare the bump shifts for different velocity thresholds for the 1x gain trials
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
velThreshLC = [0 20 40 60];
velThreshUC = [20 40 60 80];
XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for gainOneID = 1:length(allFlyData{flyID}.gainOne)
        if ~isempty(allFlyData{flyID}.gainOne{gainOneID})
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainOne{gainOneID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.gainOne{gainOneID}.RROIaveMax));
            tSpan = allFlyData{flyID}.gainOne{gainOneID}.positionDat.t - ...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(1) + ...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab)/num_planes);

            % Match the position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainOne{gainOneID}.positionDat.t >=...
                    (allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainOne{gainOneID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.gainOne{gainOneID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.gainOne{gainOneID}.positionDat.OffsetLat(tMatch(1));
            end
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(tSmooth./mean(diff(OffsetRotMatch(:,1)))));
                    
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins
                RturnsG = [];
                RturnsR = [];
                LturnsG = [];
                LturnsR = [];

                velThreshLCNow = velThreshLC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                velThreshUCNow = velThreshUC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                for i=minFG+1:maxFG
                    if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < velThreshUCNow ...
                            & abs(OffsetRotMatch(i+1-minFG,2)) < 2*pi/3
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        RturnsG = horzcat(RturnsG, GROIaveMaxFilt(:,i));
                        RturnsR = horzcat(RturnsR, RROIaveMaxFilt(:,i));
                    elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > -velThreshUCNow ...
                            & abs(OffsetRotMatch(i+1-minFG,2)) < 2*pi/3
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        LturnsG = horzcat(LturnsG, GROIaveMaxFilt(:,i));
                        LturnsR = horzcat(LturnsR, RROIaveMaxFilt(:,i));
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.gainOne{gainOneID}.RturnsG = RturnsG;
                XCVals{flyID}.vel{velCut}.gainOne{gainOneID}.RturnsR = RturnsR;
                XCVals{flyID}.vel{velCut}.gainOne{gainOneID}.LturnsG = LturnsG;
                XCVals{flyID}.vel{velCut}.gainOne{gainOneID}.LturnsR = LturnsR;
            end
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       RturnsGall = [];
       RturnsRall = [];
       LturnsGall = [];
       LturnsRall = [];
       for trial=1:length(XCVals{flyID}.vel{velCut}.gainOne)
            RturnsG = [];
            RturnsR = [];
            LturnsG = [];
            LturnsR = [];
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.gainOne{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.gainOne{trial}.RturnsG,2)
                    pkPos = find(XCVals{flyID}.vel{velCut}.gainOne{trial}.RturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.gainOne{trial}.RturnsG(:,i)));
                    RturnsG(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainOne{trial}.RturnsG(:,i),num_ROIs/2-pkPos);
                    RturnsR(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainOne{trial}.RturnsR(:,i),num_ROIs/2-pkPos);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.gainOne{trial}.LturnsG,2)
                    pkPos = find(XCVals{flyID}.vel{velCut}.gainOne{trial}.LturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.gainOne{trial}.LturnsG(:,i)));
                    LturnsG(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainOne{trial}.LturnsG(:,i),num_ROIs/2-pkPos);
                    LturnsR(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainOne{trial}.LturnsR(:,i),num_ROIs/2-pkPos);
                end
            end
            
            RturnsGall = horzcat(RturnsGall,RturnsG);
            RturnsRall = horzcat(RturnsRall,RturnsR);
            LturnsGall = horzcat(LturnsGall,LturnsG);
            LturnsRall = horzcat(LturnsRall,LturnsR);
       end
       XCVals{flyID}.vel{velCut}.RturnsGAllGainOne = RturnsGall;
       XCVals{flyID}.vel{velCut}.RturnsRAllGainOne = RturnsRall;
       XCVals{flyID}.vel{velCut}.LturnsGAllGainOne = LturnsGall;
       XCVals{flyID}.vel{velCut}.LturnsRAllGainOne = LturnsRall;
    end 
end

bumpComp = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.gainOne)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllGainOne;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllGainOne;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllGainOne;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllGainOne;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           xlim([1 16]);
           ylim([0 1.5]);
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

% set(bumpComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpComp,strcat('Results\1x_BumpComp'),'-dpdf');

%% Compare the bump shifts for different velocity thresholds for the 2x gain trials
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
velThreshLC = [0 20 40 60];
velThreshUC = [20 40 60 80];
XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo)
        if ~isempty(allFlyData{flyID}.gainTwo{gainTwoID})
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax));
            tSpan = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t - ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) + ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/num_planes);

            % Match the position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t >=...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetLat(tMatch(1));
            end
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(tSmooth./mean(diff(OffsetRotMatch(:,1)))));
                    
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins
                RturnsG = [];
                RturnsR = [];
                LturnsG = [];
                LturnsR = [];

                velThreshLCNow = velThreshLC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                velThreshUCNow = velThreshUC(velCut)*pi/180*mean(diff(OffsetRotMatch(:,1)));
                for i=minFG+1:maxFG
                    if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < velThreshUCNow ...
                            & abs(OffsetRotMatch(i+1-minFG,2)) < 2*pi/3
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        RturnsG = horzcat(RturnsG, GROIaveMaxFilt(:,i));
                        RturnsR = horzcat(RturnsR, RROIaveMaxFilt(:,i));
                    elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThreshLCNow ...
                            & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > -velThreshUCNow ...
                            & abs(OffsetRotMatch(i+1-minFG,2)) < 2*pi/3
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        LturnsG = horzcat(LturnsG, GROIaveMaxFilt(:,i));
                        LturnsR = horzcat(LturnsR, RROIaveMaxFilt(:,i));
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.gainTwo{gainTwoID}.RturnsG = RturnsG;
                XCVals{flyID}.vel{velCut}.gainTwo{gainTwoID}.RturnsR = RturnsR;
                XCVals{flyID}.vel{velCut}.gainTwo{gainTwoID}.LturnsG = LturnsG;
                XCVals{flyID}.vel{velCut}.gainTwo{gainTwoID}.LturnsR = LturnsR;
            end
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       RturnsGall = [];
       RturnsRall = [];
       LturnsGall = [];
       LturnsRall = [];
       for trial=1:length(XCVals{flyID}.vel{velCut}.gainTwo)
            RturnsG = [];
            RturnsR = [];
            LturnsG = [];
            LturnsR = [];
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.gainTwo{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.gainTwo{trial}.RturnsG,2)
                    pkPos = find(XCVals{flyID}.vel{velCut}.gainTwo{trial}.RturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.gainTwo{trial}.RturnsG(:,i)));
                    RturnsG(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainTwo{trial}.RturnsG(:,i),num_ROIs/2-pkPos);
                    RturnsR(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainTwo{trial}.RturnsR(:,i),num_ROIs/2-pkPos);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.gainTwo{trial}.LturnsG,2)
                    pkPos = find(XCVals{flyID}.vel{velCut}.gainTwo{trial}.LturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.gainTwo{trial}.LturnsG(:,i)));
                    LturnsG(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainTwo{trial}.LturnsG(:,i),num_ROIs/2-pkPos);
                    LturnsR(:,i) = circshift(XCVals{flyID}.vel{velCut}.gainTwo{trial}.LturnsR(:,i),num_ROIs/2-pkPos);
                end
            end
            
            RturnsGall = horzcat(RturnsGall,RturnsG);
            RturnsRall = horzcat(RturnsRall,RturnsR);
            LturnsGall = horzcat(LturnsGall,LturnsG);
            LturnsRall = horzcat(LturnsRall,LturnsR);
       end
       XCVals{flyID}.vel{velCut}.RturnsGAllGainTwo = RturnsGall;
       XCVals{flyID}.vel{velCut}.RturnsRAllGainTwo = RturnsRall;
       XCVals{flyID}.vel{velCut}.LturnsGAllGainTwo = LturnsGall;
       XCVals{flyID}.vel{velCut}.LturnsRAllGainTwo = LturnsRall;
    end 
end

bumpComp = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.gainTwo)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllGainTwo;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllGainTwo;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllGainTwo;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllGainTwo;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           xlim([1 16]);
           ylim([0 1.5]);
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

set(bumpComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpComp,strcat('Results\2x_BumpComp'),'-dpdf');

%% Compare PVA to gain = 2 trace
PVABehav = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)

    for gainTwoID = 1:length(allFlyData{flyID}.gainTwo)
        if ~isempty(allFlyData{flyID}.gainTwo{gainTwoID});

            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.gainTwo{gainTwoID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.gainTwo{gainTwoID}.RROIaveMax));
            tSpan = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t - ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) + ...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab)/num_planes);

            % Match the stripe position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t >=...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.OffsetLat(tMatch(1));
            end
            
            [posRot, posFor, posLat] = PositionConverter(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat);
            
            % Match the fly position data to the framegrab times
            posRotMatch = zeros(maxFG-minFG+1,1);
            posForMatch = zeros(maxFG-minFG+1,1);
            posLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t >=...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.t(1) +...
                    (allFlyData{flyID}.gainTwo{gainTwoID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
                posRotMatch(interp-minFG+1) = pi/180*posRot(tMatch(1));
                posForMatch(interp-minFG+1) = posFor(tMatch(1));
                posLatMatch(interp-minFG+1) = posLat(tMatch(1));
            end

            % Calculate the PVA for 60D05
            num_ROIs = 16;
            angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
            angsraw = angsraw';
            clear meanAngRaw;
            clear meanIntRaw;
            for ts = 1:length(GROIaveMaxFilt)
                meanAngRaw(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
                meanIntRaw(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
            end
            
            % Unwrap the values
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
            posRotUnwrap = UnWrap(posRotMatch,2,0);
            PVAUnwrap = UnWrap(meanAngRaw,1.5,1);
            
            % Plot the relative values
            subplot(length(allFlyData),7,gainTwoID+7*(flyID-1))
            hold on;
            plot(OffsetRotMatch(:,1),2*(OffsetRotUnwrap-OffsetRotUnwrap(1)),'b');
            plot(OffsetRotMatch(:,1),OffsetRotUnwrap-OffsetRotUnwrap(1),'g');
            plot(OffsetRotMatch(:,1),PVAUnwrap(minFG:maxFG)-PVAUnwrap(minFG),'k');

%             set(RGcomp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%             print(RGcomp,strcat('Results\',allFlyData{flyID}.ID,'_2x_Summary_',num2str(gainTwoID)),'-dpdf');
            
%             delete(RGcomp);
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

%% Create array of bump shapes as sorted by velocity and sign of acceleration (for dark case)
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
% velThreshLC = [0:5:90];
% velThreshUC = [5:5:95];
velThreshLC = [0:20:60];
velThreshUC = [20:20:80];
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
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(tSmooth./mean(diff(OffsetRotMatch(:,1)))));
                    
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

                for i=minFG+1:maxFG-2
                    if aFly(i+1-minFG) > 0 
                        if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < velThreshUCNow
                            RturnsGinc = horzcat(RturnsGinc, GROIaveMaxFilt(:,i));
                            RturnsRinc = horzcat(RturnsRinc, RROIaveMaxFilt(:,i));
                        elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > -velThreshUCNow
                            LturnsGinc = horzcat(LturnsGinc, GROIaveMaxFilt(:,i));
                            LturnsRinc = horzcat(LturnsRinc, RROIaveMaxFilt(:,i));
                        end
                    elseif aFly(i+1-minFG) < 0
                        if OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < velThreshUCNow
                            RturnsGdec = horzcat(RturnsGdec, GROIaveMaxFilt(:,i));
                            RturnsRdec = horzcat(RturnsRdec, RROIaveMaxFilt(:,i));
                        elseif OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) < -velThreshLCNow ...
                                & OffsetRotUnwrap(i+1-minFG) - OffsetRotUnwrap(i-minFG) > -velThreshUCNow
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

% Pool together the turns across trials
num_ROIs = 16;
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
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i)));
                    RturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(:,i)));
                    LturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(:,i),num_ROIs/2-pkPosGalign);
                    LturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGinc(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRinc(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i)));
                    RturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(:,i)));
                    LturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(:,i),num_ROIs/2-pkPosGalign);
                    LturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsGdec(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.LturnsRdec(:,i),num_ROIs/2-pkPosRalign);
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

bumpCompRinc = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc;
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
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           plot(median(LnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
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


bumpCompRdec = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec;
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
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           plot(median(LnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
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

% set(bumpCompRinc,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompRinc,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccInc'),'-dpdf');
% set(bumpCompRdec,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompRdec,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccDec'),'-dpdf');

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
% print(bumpCompRinc,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccIncEx'),'-dpdf');
print(bumpCompRinc,strcat('BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccIncEx'),'-dpdf');
set(bumpCompRdec,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompRdec,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccDecEx'),'-dpdf');
print(bumpCompRdec,strcat('BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccDecEx'),'-dpdf');

%% Get the difference between PVA direction for each fly
PVAdiff = zeros(numFlies,length(velThreshLC));
WPVAstren = zeros(numFlies,length(velThreshLC));
TPVAstren = zeros(numFlies,length(velThreshLC));
for flyID=1:numFlies
    for velCut = 1:length(velThreshLC)       
%        RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
%            XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
%            circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),-1),...
%            circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec),-1));
%        RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
%            XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
%            circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),-1),...
%            circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec),-1));
       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec));
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec));

       % Calculate the PVA differences
       num_ROIs = 16;
       angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
       angsraw = angsraw';
       PVAdiff(flyID,velCut) = -circ_mean(angsraw,mean(RturnsR,2))+circ_mean(angsraw,mean(RturnsG,2));
       WPVAstren(flyID,velCut) = circ_r(angsraw,mean(RturnsR,2));
       TPVAstren(flyID,velCut) = circ_r(angsraw,mean(RturnsG,2));
       PVAdiff(flyID,velCut) = min(PVAdiff(flyID,velCut),2*pi-PVAdiff(flyID,velCut));
       
    end

end

PVAdiff = -PVAdiff;

PVAThresh = 0.12;

PVAdiffplot = figure;
set(gcf,'Position',[100 0 560 420*2.5]);
subplot(6,1,[1:4]);
hold on;
for PVAplt = 1:6
    PVAcolor = zeros(numFlies,3);
    for clrNow = 1:numFlies
        cNow = 1-TPVAstren(clrNow,PVAplt)/PVAThresh;
        PVAcolor(clrNow,:) = [cNow cNow cNow];
    end
    scatter(zeros(numFlies,1)+PVAplt-1,180/pi*PVAdiff(:,PVAplt),50,PVAcolor,'filled');
    line([-1.5+PVAplt -0.5+PVAplt],[180/pi*mean(PVAdiff(:,PVAplt)) 180/pi*mean(PVAdiff(:,PVAplt))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',10);
ylabel('PVA difference (deg)');
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('vRot bins (o/s)');
xlim([-0.5 5.5]);
ylim([-10 40]);

% hold on;
% scatter(zeros(numFlies,1),180/pi*PVAdiff(:,1),'filled','k');
% line([-0.5 0.5],[180/pi*mean(PVAdiff(:,1)) 180/pi*mean(PVAdiff(:,1))],'Color','k','LineWidth',2);
% scatter(zeros(numFlies,1)+1,180/pi*PVAdiff(:,2),'filled','k');
% line([0.5 1.5],[180/pi*mean(PVAdiff(:,2)) 180/pi*mean(PVAdiff(:,2))],'Color','k','LineWidth',2);
% scatter(zeros(numFlies,1)+2,180/pi*PVAdiff(:,3),'filled','k');
% line([1.5 2.5],[180/pi*mean(PVAdiff(:,3)) 180/pi*mean(PVAdiff(:,3))],'Color','k','LineWidth',2);
% scatter(zeros(numFlies,1)+3,180/pi*PVAdiff(:,4),'filled','k');
% line([2.5 3.5],[180/pi*mean(PVAdiff(:,4)) 180/pi*mean(PVAdiff(:,4))],'Color','k','LineWidth',2);
% scatter(zeros(numFlies,1)+4,180/pi*PVAdiff(:,5),'filled','k');
% line([3.5 4.5],[180/pi*mean(PVAdiff(:,5)) 180/pi*mean(PVAdiff(:,5))],'Color','k','LineWidth',2);
% scatter(zeros(numFlies,1)+5,180/pi*PVAdiff(:,6),'filled','k');
% line([4.5 5.5],[180/pi*mean(PVAdiff(:,6)) 180/pi*mean(PVAdiff(:,6))],'Color','k','LineWidth',2);
% set(gca,'FontSize',14);
% ylabel('PVA difference (deg)');
% set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
% xlabel('v_r bins (deg/sec)');
% xlim([-0.5 5.5]);
% ylim([-10 40]);
% ylim([-5 30]);
subplot(6,1,5);
imagesc(mean(TPVAstren,1));
caxis([0 0.12]);
colorbar;
colormap(brewermap(64, 'Greys'));
xlim([0.5 6.5]);

subplot(6,1,6);
imagesc(mean(WPVAstren,1));
caxis([0 0.12]);
colorbar;
colormap(brewermap(64, 'Greys'));
xlim([0.5 6.5]);

set(PVAdiffplot,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(PVAdiffplot,strcat('Results\','PVADiffnStrenOtherColorDir'),'-dpdf');
    
% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,2),PVAdiff(:,3)),vertcat(zeros(10,1), ones(10,1), 2*ones(10,1)))
% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,2)),vertcat(zeros(5,1), ones(5,1)))
% anova1(vertcat(PVAdiff(:,2),PVAdiff(:,3)),vertcat(zeros(5,1), ones(5,1)))
% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,3)),vertcat(zeros(5,1), ones(5,1)))

%% Get the peak heights and widths 
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=[1:3]
    for velCut = 1:length(velThreshLC)     
       GturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec));
       GturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec));
       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec));
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec));

       if ~isempty(RturnsG) 
           % Calculate the PVA differences
           GpkMax(flyID, velCut) = max(mean(RturnsG,2))-min(mean(RturnsG,2));
           GpkRng = find(mean(RturnsG,2)-min(mean(RturnsG,2))>0.5*GpkMax(flyID, velCut));
           if ~isempty(GpkRng)
               GpkWid(flyID, velCut) = GpkRng(end)-GpkRng(1);
           end

           RpkMax(flyID, velCut) = max(mean(GturnsR,2))-min(mean(GturnsR,2));
           RpkRng = find(mean(GturnsR,2)-min(mean(GturnsR,2))>0.5*RpkMax(flyID, velCut));
           if ~isempty(RpkRng)
               RpkWid(flyID, velCut) = RpkRng(end)-RpkRng(1);       
           end
       end
    end

end


GNorm = GpkMax(:, 1);
RNorm = RpkMax(:, 1);
for flyID = 1:size(GpkMax,1)
    for velCut = 1:length(velThreshLC)
        GpkMax(flyID,velCut) = GpkMax(flyID,velCut)/GNorm(flyID);
        RpkMax(flyID,velCut) = RpkMax(flyID,velCut)/RNorm(flyID);
    end
end
        

PVAdiff = PVAdiff;

figure;
subplot(2,2,1);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+30*i-15,squeeze(GpkMax(:,i)),'filled','k'); 
   line([30*i-30 30*i],[mean(GpkMax(:,i)) mean(GpkMax(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 185]);
ylim([0.5 8]);
xlabel('vr bins (deg/sec)');

subplot(2,2,2);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+30*i-15,22.5*squeeze(GpkWid(:,i)),'filled','k'); 
   line([30*i-30 30*i],[22.5*mean(GpkWid(:,i)) 22.5*mean(GpkWid(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 185]);
xlabel('vr bins (deg/sec)');
ylim([0 250]);

subplot(2,2,3);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+30*i-15,squeeze(RpkMax(:,i)),'filled','k'); 
   line([30*i-30 30*i],[mean(RpkMax(:,i)) mean(RpkMax(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 185]);
xlabel('vr bins (deg/sec)');

subplot(2,2,4);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+30*i-15,22.5*squeeze(RpkWid(:,i)),'filled','k'); 
   line([30*i-30 30*i],[22.5*mean(RpkWid(:,i)) 22.5*mean(RpkWid(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 185]);
xlabel('vr bins (deg/sec)');
ylim([0 250]);


%% Align the peaks of 37F06, align to one peak, and plot on a circle
trial = 'dark';

if strcmp(redLine,'37F06')
    tileDataName = 'RROIaveMax';
    wedgeDataName = 'GROIaveMax';
else
    tileDataName = 'GROIaveMax';
    wedgeDataName = 'RROIaveMax';
end
    
% Specify the range of rotational velocities to consider
vRMin = 0*pi/180;
vRMax = 360*pi/180;
vRSpan = 30*pi/180;
vRBinNum = round((vRMax-vRMin)/vRSpan);

tileData = cell(size(allFlyData,2));
wedgeData = cell(size(allFlyData,2));
% Sort the data
for flyID = 1:length(allFlyData)
    tileData{flyID}.CW = cell(vRBinNum,1);
    tileData{flyID}.CCW = cell(vRBinNum,1);
    tileData{flyID}.Stop = [];
    
    for trialID = 1:length(allFlyData{flyID}.(trial))
        vRot = allFlyData{flyID}.(trial){trialID}.positionDatMatch.vRot;
        Sig = allFlyData{flyID}.(trial){trialID}.(tileDataName);
        for tStep = 1:length(vRot)
            pkMax = max(Sig(:,tStep+1));
            pkShift = -find(Sig(:,tStep+1) == pkMax);
            SigNow = circshift(Sig(:,tStep+1),pkShift);
            
            if vRot(tStep) == 0
                tileData{flyID}.Stop = horzcat(tileData{flyID}.Stop, SigNow);
            elseif vRot(tStep) > 0
                tileData{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(tileData{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    SigNow);
            else
                tileData{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(tileData{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    SigNow);
            end
        end
    end
    
    wedgeData{flyID}.CW = cell(vRBinNum,1);
    wedgeData{flyID}.CCW = cell(vRBinNum,1);
    wedgeData{flyID}.Stop = [];
    
    for trialID = 1:length(allFlyData{flyID}.(trial))
        vRot = allFlyData{flyID}.(trial){trialID}.positionDatMatch.vRot;
        Sig = allFlyData{flyID}.(trial){trialID}.(wedgeDataName);
        for tStep = 1:length(vRot)
            pkMax = max(Sig(:,tStep+1));
            pkShift = -find(Sig(:,tStep+1) == pkMax);
            SigNow = circshift(Sig(:,tStep+1),pkShift);
            
            if vRot(tStep) == 0
                wedgeData{flyID}.Stop = horzcat(wedgeData{flyID}.Stop, SigNow);
            elseif vRot(tStep) > 0
                wedgeData{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(wedgeData{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    SigNow);
            else
                wedgeData{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(wedgeData{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    SigNow);
            end
        end
    end
end


circCent = [128 128];
circHeight = 100;
circWidth = 100;
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;

for flyID = 1:length(allFlyData)
    bumpWidthComp = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(4,7,1);
    title(strcat('tile EB activity, fly ',num2str(flyID)));
    hold on;
    stopData = mean(tileData{flyID}.Stop,2);
    stopData = (stopData-min(stopData))./(max(stopData)-min(stopData));
    for incROI = 1:num_ROIs
        xWedge(1) = circCent(1);
        yWedge(1) = circCent(2);
        for angs=1:5
            xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        end
        patch(xWedge,yWedge,stopData(incROI));
    end
    axis off;
    axis equal;
            
    for bin=1:6
        subplot(4,7,1+bin);
        title(strcat('CW, ', num2str(vRMin+vRSpan*(bin-1)),'-',num2str(vRMin+vRSpan*bin),'rad/sec'));
        CWData = mean(tileData{flyID}.CW{bin},2);
        CWData = (CWData-min(CWData))./(max(CWData)-min(CWData));
        if ~isempty(CWData)
            for incROI = 1:num_ROIs
               xWedge(1) = circCent(1);
            yWedge(1) = circCent(2);
            for angs=1:5
                xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            end
            patch(xWedge,yWedge,CWData(incROI));
            end
            axis off;
            axis equal;
        end
        
        subplot(4,7,8+bin);
        title(strcat('CCW, ', num2str(vRMin+vRSpan*(bin-1)),'-',num2str(vRMin+vRSpan*bin),'rad/sec'));
        CCWData = mean(tileData{flyID}.CCW{bin},2);
        CCWData = (CCWData-min(CCWData))./(max(CCWData)-min(CCWData));
        if ~isempty(CCWData)
            for incROI = 1:num_ROIs
               xWedge(1) = circCent(1);
            yWedge(1) = circCent(2);
            for angs=1:5
                xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            end
            patch(xWedge,yWedge,CCWData(incROI));
            end
            axis off;
            axis equal;
        end
    end
    
    subplot(4,7,15);
    title(strcat('wedge EB activity, fly ',num2str(flyID)));
    hold on;
    stopData = mean(wedgeData{flyID}.Stop,2);
    stopData = (stopData-min(stopData))./(max(stopData)-min(stopData));
    for incROI = 1:num_ROIs
        xWedge(1) = circCent(1);
        yWedge(1) = circCent(2);
        for angs=1:5
            xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
        end
        patch(xWedge,yWedge,stopData(incROI));
    end
    axis off;
    axis equal;
            
    for bin=1:6
        subplot(4,7,15+bin);
        title(strcat('CW, ', num2str(vRMin+vRSpan*(bin-1)),'-',num2str(vRMin+vRSpan*bin),'rad/sec'));
        CWData = mean(wedgeData{flyID}.CW{bin},2);
        CWData = (CWData-min(CWData))./(max(CWData)-min(CWData));
        if ~isempty(CWData)
            for incROI = 1:num_ROIs
               xWedge(1) = circCent(1);
            yWedge(1) = circCent(2);
            for angs=1:5
                xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            end
            patch(xWedge,yWedge,CWData(incROI));
            end
            axis off;
            axis equal;
        end
        
        subplot(4,7,22+bin);
        title(strcat('CCW, ', num2str(vRMin+vRSpan*(bin-1)),'-',num2str(vRMin+vRSpan*bin),'rad/sec'));
        CCWData = mean(wedgeData{flyID}.CCW{bin},2);
        CCWData = (CCWData-min(CCWData))./(max(CCWData)-min(CCWData));
        if ~isempty(CCWData)
            for incROI = 1:num_ROIs
               xWedge(1) = circCent(1);
            yWedge(1) = circCent(2);
            for angs=1:5
                xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            end
            patch(xWedge,yWedge,CCWData(incROI));
            end
            axis off;
            axis equal;
        end
    end
    colormap(brewermap(64, 'Blues'));
    set(bumpWidthComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(bumpWidthComp,strcat('Results\','Fly',num2str(flyID),'_NormImBumpComp_',redLine,'Red_',greenLine,'Green'),'-dpdf');
end

%% Create array of bump shapes as sorted by velocity and sign of acceleration - fixed for new data analysis

visDisp = 'dark';
flyRange = 1:4;%length(allFlyData);
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
for flyID = flyRange
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.(visDisp))
        if ~isempty(allFlyData{flyID}.(visDisp){darkID})
                 
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins by increasing or
                % decreasing velocity
                RturnsG = [];
                RturnsR = [];
                LturnsG = [];
                LturnsR = [];

                velThreshUCNow = velThreshUC(velCut)*pi/180;
                velThreshLCNow = velThreshLC(velCut)*pi/180;

                minFG = allFlyData{flyID}.(visDisp){darkID}.positionDatMatch.minFG;
                maxFG = allFlyData{flyID}.(visDisp){darkID}.positionDatMatch.maxFG;
                vRot = allFlyData{flyID}.(visDisp){darkID}.positionDatMatch.vRot;
                
                for i=1:maxFG-minFG
                    if vRot(i) > velThreshLCNow ...
                            & vRot(i) < velThreshUCNow
                        RturnsG = horzcat(RturnsG, allFlyData{flyID}.(visDisp){darkID}.GROIaveMax(:,i));
                        RturnsR = horzcat(RturnsR, allFlyData{flyID}.(visDisp){darkID}.RROIaveMax(:,i));
                    end
                    if -vRot(i) > velThreshLCNow ...
                            & -vRot(i) < velThreshUCNow
                        LturnsG = horzcat(LturnsG, allFlyData{flyID}.(visDisp){darkID}.GROIaveMax(:,i));
                        LturnsR = horzcat(LturnsR, allFlyData{flyID}.(visDisp){darkID}.RROIaveMax(:,i));
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.(visDisp){darkID}.RturnsG = RturnsG;
                XCVals{flyID}.vel{velCut}.(visDisp){darkID}.RturnsR = RturnsR;
                XCVals{flyID}.vel{velCut}.(visDisp){darkID}.LturnsG = LturnsG;
                XCVals{flyID}.vel{velCut}.(visDisp){darkID}.LturnsR = LturnsR;
            end
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = flyRange
    for velCut = 1:length(velThreshLC)       
       RturnsGallGalign = [];
       RturnsRallGalign = [];
       LturnsGallGalign = [];
       LturnsRallGalign = [];
       
       RturnsGallRalign = [];
       RturnsRallRalign = [];
       LturnsGallRalign = [];
       LturnsRallRalign = [];
        for trial=1:length(XCVals{flyID}.vel{velCut}.(visDisp))
            RturnsGGalign = [];
            RturnsRGalign = [];
            LturnsGGalign = [];
            LturnsRGalign = [];
            
            RturnsGRalign = [];
            RturnsRRalign = [];
            LturnsGRalign = [];
            LturnsRRalign = [];
            
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.(visDisp){trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsG(:,i)));
                    RturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsG(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsR(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsG(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.RturnsR(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsG(:,i)));
                    LturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsG(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsR(:,i),num_ROIs/2-pkPosGalign);
                    LturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsG(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visDisp){trial}.LturnsR(:,i),num_ROIs/2-pkPosRalign);
                end
            end
            
            RturnsGallGalign = horzcat(RturnsGallGalign,RturnsGGalign);
            RturnsRallGalign = horzcat(RturnsRallGalign,RturnsRGalign);
            LturnsGallGalign = horzcat(LturnsGallGalign,LturnsGGalign);
            LturnsRallGalign = horzcat(LturnsRallGalign,LturnsRGalign);
            
            RturnsGallRalign = horzcat(RturnsGallRalign,RturnsGRalign);
            RturnsRallRalign = horzcat(RturnsRallRalign,RturnsRRalign);
            LturnsGallRalign = horzcat(LturnsGallRalign,LturnsGRalign);
            LturnsRallRalign = horzcat(LturnsRallRalign,LturnsRRalign);
            
       end
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalign = RturnsGallGalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalign = RturnsRallGalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalign = LturnsGallGalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalign = LturnsRallGalign;
       
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalign = RturnsGallRalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalign = RturnsRallRalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalign = LturnsGallRalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalign = LturnsRallRalign;
    end 
end

bumpComp = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = flyRange
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.(visDisp))
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalign;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           
           subplot(2*length(flyRange),length(velThreshLC),2*(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           plot(quantile(LnetR,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetR,0.75,2),'m','LineWidth',1);
           plot(median(LnetR,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetR,2)),'Color','m');
           plot(median(LnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
           set(gca,'FontSize',14);
           
           if flyID == 1
               title(strcat('v_{thresh} = ',num2str(XCVals{flyID}.vel{velCut}.threshold),'deg/sec'));
           end
           
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalign;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
           LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
          
           subplot(2*length(flyRange),length(velThreshLC),(2*flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetG,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetG,0.75,2),'b','LineWidth',1);
           plot(median(RnetG,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetG,2)),'Color','b');
           plot(median(RnetR,2),'--r','LineWidth',1);
           
           plot(quantile(LnetG,0.25,2),'m','LineWidth',1);
           plot(quantile(LnetG,0.75,2),'m','LineWidth',1);
           plot(median(LnetG,2),'m','LineWidth',4);
           text(15,0.8,num2str(size(LnetG,2)),'Color','m');
           plot(median(LnetR,2),'--r','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
           set(gca,'FontSize',14);
       end
       if flyID == numFlies
           xlabel('ROI');
       end
    end
    subplot(2*length(flyRange),length(velThreshLC),2*(flyID-1)*length(velThreshLC)+1);
    ylabel(allFlyData{flyID}.ID);
end

%% Create array of bump shapes as sorted by velocity and sign of acceleration (for dark case)
visStim = 'dark';

% Threshold velocities to consider
velThreshLC = [0:30:150];
velThreshUC = [30:30:180];
% velThreshLC = [0:45:135];
% velThreshUC = [45:45:180];
% velThreshLC = [0:60:120];
% velThreshUC = [60:60:180];
% velThreshLC = [-720];
% velThreshUC = [720];
% velThreshLC = [0];
% velThreshUC = [10];

XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.(visStim))
        if ~isempty(allFlyData{flyID}.(visStim){darkID}) & ~(flyID ==4 & darkID ==1)
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
                
                
                minFG = allFlyData{flyID}.(visStim){darkID}.positionDatMatch.minFG;
                maxFG = allFlyData{flyID}.(visStim){darkID}.positionDatMatch.maxFG;
                GROIaveMax = allFlyData{flyID}.(visStim){darkID}.GROIaveMax;
                RROIaveMax = allFlyData{flyID}.(visStim){darkID}.RROIaveMax;
                vR = allFlyData{flyID}.(visStim){darkID}.positionDatMatch.vRot*180/pi;
                aFly = diff(vR);

                for i=1:length(GROIaveMax)-2
                    if aFly(i) > 0 
                        if vR(i) > velThreshLC(velCut) ...
                                & vR(i) < velThreshUC(velCut)
                            RturnsGinc = horzcat(RturnsGinc, GROIaveMax(:,i));
                            RturnsRinc = horzcat(RturnsRinc, RROIaveMax(:,i));
                        elseif vR(i) < -velThreshLC(velCut) ...
                                & vR(i) > -velThreshUC(velCut)
                            LturnsGinc = horzcat(LturnsGinc, GROIaveMax(:,i));
                            LturnsRinc = horzcat(LturnsRinc, RROIaveMax(:,i));
                        end
                    elseif aFly(i) < 0
                        if vR(i) > velThreshLC(velCut) ...
                                & vR(i) < velThreshUC(velCut)
                            RturnsGdec = horzcat(RturnsGdec, GROIaveMax(:,i));
                            RturnsRdec = horzcat(RturnsRdec, RROIaveMax(:,i));
                        elseif vR(i) < -velThreshLC(velCut) ...
                                & vR(i) > -velThreshUC(velCut)
                            LturnsGdec = horzcat(LturnsGdec, GROIaveMax(:,i));
                            LturnsRdec = horzcat(LturnsRdec, RROIaveMax(:,i));
                        end
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.RturnsGinc = RturnsGinc;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.RturnsGdec = RturnsGdec;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.RturnsRinc = RturnsRinc;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.RturnsRdec = RturnsRdec;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.LturnsGinc = LturnsGinc;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.LturnsGdec = LturnsGdec;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.LturnsRinc = LturnsRinc;
                XCVals{flyID}.vel{velCut}.(visStim){darkID}.LturnsRdec = LturnsRdec;
            end
            
            clear RROIaveMax GROIaveMax vR;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
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
        for trial=1:length(XCVals{flyID}.vel{velCut}.(visStim))
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
            if ~isempty(XCVals{flyID}.vel{velCut}.(visStim){trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRinc(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGinc(:,i)));
                    RturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGinc(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRinc(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRinc(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGinc(:,i)));
                    LturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGinc(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRinc(:,i),num_ROIs/2-pkPosGalign);
                    LturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGinc(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRinc(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRdec(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGdec(:,i)));
                    RturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsGdec(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.RturnsRdec(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRdec(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGdec(:,i)));
                    LturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGdec(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRdec(:,i),num_ROIs/2-pkPosGalign);
                    LturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsGdec(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.(visStim){trial}.LturnsRdec(:,i),num_ROIs/2-pkPosRalign);
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

%% Plot the above for a given fly, pulling out the two directions separately.
flyID = 7;
visStim = 'dark';
bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);

for velCut = 1:length(velThreshLC)       
   for trial=1:length(XCVals{flyID}.vel{velCut}.(visStim))
       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;

       subplot(8,length(velThreshLC),velCut);
       hold on;
       imagesc(mean(RnetG,2)');
       axis off;
       caxis([0 1.5]);
       colorbar;
       title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));

       subplot(8,length(velThreshLC),length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(RnetR,2)');
       axis off;
       caxis([0 2]);
       colorbar;
       title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));
       
       subplot(8,length(velThreshLC),2*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(LnetG,2)');
       axis off;
       caxis([0 1.5]);
       colorbar;
       
       subplot(8,length(velThreshLC),3*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(LnetR,2)');
       axis off;
       caxis([0 2]);
       colorbar;

       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;


       subplot(8,length(velThreshLC),4*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(RnetG,2)');
       axis off;
       caxis([0 1.5]);
       colorbar;

       subplot(8,length(velThreshLC),5*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(RnetR,2)');
       axis off;
       caxis([0 2]);
       colorbar;
       
       subplot(8,length(velThreshLC),6*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(LnetG,2)');
       axis off;
       caxis([0 1.5]);
       colorbar;
       
       subplot(8,length(velThreshLC),7*length(velThreshLC)+velCut);
       hold on;
       imagesc(mean(LnetR,2)');
       axis off;
       caxis([0 2]);
       colorbar;

       
   end
end

subplot(8,length(velThreshLC),1);
title('acc > 0');
axis off;

subplot(8,length(velThreshLC),3*length(velThreshLC)+1);
title('acc < 0');
axis off;

% colormap(brewermap(64, 'Blues'));
blues = colormap(brewermap(64, 'Blues'));
greens = blues;
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
colormap(greens);
% reds = blues;
% reds(:,1) = blues(:,3);
% reds(:,2) = blues(:,1);
% reds(:,3) = blues(:,1);
% colormap(reds);

% set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompR,strcat('BumpCompGAlignNormPlots_Red',redLine,'_Green',greenLine,'AccIncEx_2_greens'),'-dpdf');

%% Plot the above for a given fly, pulling out the two directions separately, using circular plotting.
flyID = 7;
visStim = 'dark';
bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);

circCent = [128 128];
circHeight = 100;
circWidth = 100;
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;

num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs;
angsraw = angsraw';

for velCut = 1:length(velThreshLC)       
   for trial=1:length(XCVals{flyID}.vel{velCut}.(visStim))
       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;

       subplot(8,length(velThreshLC),velCut);
       hold on;
       RGPlot = mean(RnetG,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RGPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 1]);
       colorbar;
       title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));
       PVAang = circ_mean(angsraw, RGPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','g');

       subplot(8,length(velThreshLC),length(velThreshLC)+velCut);
       hold on;
       RRPlot = mean(RnetR,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RRPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 2]);
       colorbar;
       PVAang = circ_mean(angsraw, RRPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','r');

       
       subplot(8,length(velThreshLC),2*length(velThreshLC)+velCut);
       hold on;
       LGPlot = mean(LnetG,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LGPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 1]);
       colorbar;
       PVAang = circ_mean(angsraw, LGPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','g');
       
       subplot(8,length(velThreshLC),3*length(velThreshLC)+velCut);
       hold on;
       LRPlot = mean(LnetR,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LRPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 2]);
       colorbar;
       PVAang = circ_mean(angsraw, LRPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','r');

       RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec;
       RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec;
       LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec;
       LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec;
       RnetR = RturnsR-1;
       RnetG = RturnsG-1;
       LnetR = LturnsR-1;
       LnetG = LturnsG-1;


       subplot(8,length(velThreshLC),4*length(velThreshLC)+velCut);
       hold on;
       RGPlot = mean(RnetG,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RGPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 1]);
       colorbar;
       title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));
       PVAang = circ_mean(angsraw, RGPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','g');

       subplot(8,length(velThreshLC),5*length(velThreshLC)+velCut);
       hold on;
       RRPlot = mean(RnetR,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RRPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 2]);
       colorbar;
       PVAang = circ_mean(angsraw, RRPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','r');
       
       subplot(8,length(velThreshLC),6*length(velThreshLC)+velCut);
       hold on;
       LGPlot = mean(LnetG,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LGPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 1]);
       colorbar;
       PVAang = circ_mean(angsraw, LGPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','g');
       
       subplot(8,length(velThreshLC),7*length(velThreshLC)+velCut);
       hold on;
       LRPlot = mean(LnetR,2);
       for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LRPlot(incROI));
       end
       axis off;
       axis equal;
       caxis([0 2]);
       colorbar;
       PVAang = circ_mean(angsraw, LRPlot);
       line([circCent(1) circCent(1)+circWidth*cos(PVAang)], [circCent(2) circCent(2)+circHeight*sin(PVAang)],'Color','r');

       
   end
end

subplot(8,length(velThreshLC),1);
title('acc > 0');
axis off;

subplot(8,length(velThreshLC),3*length(velThreshLC)+1);
title('acc < 0');
axis off;

% colormap(brewermap(64, 'Blues'));
blues = colormap(brewermap(64, 'Blues'));
greens = blues;
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
colormap(greens);
% reds = blues;
% reds(:,1) = blues(:,3);
% reds(:,2) = blues(:,1);
% reds(:,3) = blues(:,1);
% colormap(reds);

set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompR,strcat('BumpCompGAlignNormPlots_Red',redLine,'_Green',greenLine,'AccIncEx_4_greens'),'-dpdf');

%% Plot the above for a given fly, pulling out the two directions separately, using circular plotting, condensed.
flyID = 10;
visStim = 'dark';
bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);

circCent = [128 128];
circHeight = 100;
circWidth = 100;
num_ROIs = 16;
deltaAng = 2*pi/num_ROIs;

num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs;
angsraw = angsraw';

for velCut = 1:length(velThreshLC)
    if velCut ~= length(velThreshLC)
        RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc;
        RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc;
        LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc;
        LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc;
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;

        subplot(4,2*length(velThreshLC)-1,velCut);
        hold on;
        RGPlot = mean(RnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));

        subplot(4,2*length(velThreshLC)-1,2*length(velThreshLC)-1+velCut);
        hold on;
        RRPlot = mean(RnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, RGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, RRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');


        subplot(4,2*length(velThreshLC)-1,4*length(velThreshLC)-2+velCut);
        hold on;
        LGPlot = mean(LnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;

        subplot(4,2*length(velThreshLC)-1,6*length(velThreshLC)-3+velCut);
        hold on;
        LRPlot = mean(LnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, LGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, LRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');

        RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec;
        RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec;
        LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec;
        LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec;
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;


        subplot(4,2*length(velThreshLC)-1,2*length(velThreshLC)-1-velCut+1);
        hold on;
        RGPlot = mean(RnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));

        subplot(4,2*length(velThreshLC)-1,4*length(velThreshLC)-2-velCut+1);
        hold on;
        RRPlot = mean(RnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, RGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, RRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');

        subplot(4,2*length(velThreshLC)-1,6*length(velThreshLC)-3-velCut+1);
        hold on;
        LGPlot = mean(LnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;

        subplot(4,2*length(velThreshLC)-1,8*length(velThreshLC)-4-velCut+1);
        hold on;
        LRPlot = mean(LnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, LGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, LRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');
    else
        RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec);
        RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec);
        LturnsG = horzcat(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec);
        LturnsR = horzcat(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec);
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;

        subplot(4,2*length(velThreshLC)-1,velCut);
        hold on;
        RGPlot = mean(RnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));

        subplot(4,2*length(velThreshLC)-1,2*length(velThreshLC)-1+velCut);
        hold on;
        RRPlot = mean(RnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,RRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, RGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, RRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');


        subplot(4,2*length(velThreshLC)-1,4*length(velThreshLC)-2+velCut);
        hold on;
        LGPlot = mean(LnetG,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LGPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;

        subplot(4,2*length(velThreshLC)-1,6*length(velThreshLC)-3+velCut);
        hold on;
        LRPlot = mean(LnetR,2);
        for incROI = 1:num_ROIs
           xWedge(1) = circCent(1);
           yWedge(1) = circCent(2);
           for angs=1:5
               xWedge(angs+1) = circCent(1)+circWidth/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
               yWedge(angs+1) = circCent(2)+circHeight/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
           end
           patch(xWedge,yWedge,LRPlot(incROI));
        end
        axis off;
        axis equal;
        caxis([0 1]);
        colorbar;
        PVAang = circ_mean(angsraw, LGPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','g');
        PVAang = circ_mean(angsraw, LRPlot)-deltaAng/2;
        line([circCent(1) circCent(1)+circWidth/2*cos(PVAang)], [circCent(2) circCent(2)+circHeight/2*sin(PVAang)],'Color','r');
        
    end
end

% colormap(brewermap(64, 'Blues'));
blues = colormap(brewermap(64, 'Blues'));
% greens = blues;
% greens(:,2) = blues(:,3);
% greens(:,3) = blues(:,2);
% colormap(greens);
% reds = blues;
% reds(:,1) = blues(:,3);
% reds(:,2) = blues(:,1);
% reds(:,3) = blues(:,1);
% colormap(reds);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
colormap(magentas);


set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompR,strcat('BumpCompGAlignNormPlots_Red',redLine,'_Green',greenLine,'AccIncEx_6_magentas'),'-dpdf');

%% Plot the above for a given fly, pulling out the two directions separately, using cross sections.
flyID = 10;
visStim = 'dark';
bumpCompR = figure('units','normalized','outerposition',[0 0 1 1]);


num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs;
angsraw = angsraw';

for velCut = 1:length(velThreshLC)
    if velCut ~= length(velThreshLC)
        RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc;
        RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc;
        LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc;
        LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc;
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;

        subplot(2,2*length(velThreshLC)-1,velCut);
        hold on;
        plot(angsraw,mean(RnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(RnetG,2)-std(RnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(RnetG,2)+std(RnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(RnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(RnetR,2)-std(RnetG,[],2),'m','LineWidth',1);
        plot(angsraw,mean(RnetR,2)+std(RnetG,[],2),'m','LineWidth',1);
        
        PVAang = 2*pi-abs(circ_mean(angsraw, mean(RnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = abs(circ_mean(angsraw, mean(RnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');
        
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));

        subplot(2,2*length(velThreshLC)-1,2*length(velThreshLC)-1+velCut);
        hold on;
        plot(angsraw,mean(LnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(LnetG,2)-std(LnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(LnetG,2)+std(LnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(LnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(LnetR,2)-std(LnetG,[],2),'m','LineWidth',1);
        plot(angsraw,mean(LnetR,2)+std(LnetG,[],2),'m','LineWidth',1);
        
        PVAang = abs(circ_mean(angsraw, mean(LnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = 2*pi-abs(circ_mean(angsraw, mean(LnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');

        RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec;
        RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec;
        LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec;
        LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec;
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');

        subplot(2,2*length(velThreshLC)-1,2*length(velThreshLC)-velCut);
        hold on;
        plot(angsraw,mean(RnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(RnetG,2)-std(RnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(RnetG,2)+std(RnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(RnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(RnetR,2)-std(RnetR,[],2),'m','LineWidth',1);
        plot(angsraw,mean(RnetR,2)+std(RnetR,[],2),'m','LineWidth',1);
        
        PVAang = 2*pi-abs(circ_mean(angsraw, mean(RnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = abs(circ_mean(angsraw, mean(RnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');
        
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');

        subplot(2,2*length(velThreshLC)-1,4*length(velThreshLC)-1-velCut);
        hold on;
        plot(angsraw,mean(LnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(LnetG,2)-std(LnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(LnetG,2)+std(LnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(LnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(LnetR,2)-std(LnetR,[],2),'m','LineWidth',1);
        plot(angsraw,mean(LnetR,2)+std(LnetR,[],2),'m','LineWidth',1);
        
        PVAang = abs(circ_mean(angsraw, mean(LnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = abs(circ_mean(angsraw, mean(LnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');
        
    else
        RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec);
        RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec);
        LturnsG = horzcat(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec);
        LturnsR = horzcat(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc,XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec);
        RnetR = RturnsR-1;
        RnetG = RturnsG-1;
        LnetR = LturnsR-1;
        LnetG = LturnsG-1;

        subplot(2,2*length(velThreshLC)-1,2*length(velThreshLC)-velCut);
        hold on;
        plot(angsraw,mean(RnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(RnetG,2)-std(RnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(RnetG,2)+std(RnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(RnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(RnetR,2)-std(RnetR,[],2),'m','LineWidth',1);
        plot(angsraw,mean(RnetR,2)+std(RnetR,[],2),'m','LineWidth',1);
        
        PVAang = 2*pi-abs(circ_mean(angsraw, mean(RnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = abs(circ_mean(angsraw, mean(RnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');
        
        title(strcat('v_{rot} = ',num2str(velThreshLC(velCut)),'-',num2str(velThreshUC(velCut)),'deg/sec'));
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');

        subplot(2,2*length(velThreshLC)-1,4*length(velThreshLC)-1-velCut);
        hold on;
        plot(angsraw,mean(LnetG,2),'g','LineWidth',2);
        plot(angsraw,mean(LnetG,2)-std(LnetG,[],2),'g','LineWidth',1);
        plot(angsraw,mean(LnetG,2)+std(LnetG,[],2),'g','LineWidth',1);

        plot(angsraw,mean(LnetR,2),'m','LineWidth',2);
        plot(angsraw,mean(LnetR,2)-std(LnetR,[],2),'m','LineWidth',1);
        plot(angsraw,mean(LnetR,2)+std(LnetR,[],2),'m','LineWidth',1);
        
        PVAang = abs(circ_mean(angsraw, mean(LnetG,2)));
        line([PVAang PVAang],[0 1],'Color','g');
        PVAang = abs(circ_mean(angsraw, mean(LnetR,2)));
        line([PVAang PVAang],[0 1],'Color','m');
        
        xlim([0 2*pi]);
        ylim([-0.25 1.5]);
        xlabel('EB position (rad)');
        ylabel('\DeltaF/F');
    end
end

set(bumpCompR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpCompR,strcat('BumpCompGAlignNormPlots_Red',redLine,'_Green',greenLine,'AccIncEx_Profile'),'-dpdf');

%% Get and plot the correlation coefficients between the PVA and the fly's rotational position for each fly
trial = 'dark';
tCut = 15;
CCs = zeros(length(allFlyData),5);
for flyID = 1:length(allFlyData)
    for darkID = 1:length(allFlyData{flyID}.(trial))
        % Get the unwrapped rotational position
        Rot = UnWrap(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,2),2,0);
        sgolayOrder = 3;
        sgolayWindow = 11;
        Rot = sgolayfilt(Rot,sgolayOrder,sgolayWindow);
            
        % Find the unwrapped PVA
        RSig = allFlyData{flyID}.(trial){darkID}.RROIaveMax;
        
        % Calculate the PVA for 60D05
        num_ROIs = 16;
        angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
        angsraw = angsraw';
        clear meanAngRaw;
        clear meanIntRaw;
        for ts = 1:length(RSig)
            meanAngRaw(ts) = circ_mean(angsraw, squeeze(RSig(:,ts)));
            meanIntRaw(ts) = circ_r(angsraw, squeeze(RSig(:,ts)));
        end
        
        PVA = UnWrap(meanAngRaw,1.5,0);
        
        tRange = find(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,1) > tCut);
        CC = corrcoef(PVA(tRange),Rot(tRange));
        CCs(flyID,darkID) = CC(1,2);
        
%         if CC(1,2) < 0
%             figure;
%             hold on;
%             plot(Rot);
%             plot(PVA);
%         end
    end
end
    
figure;
hold on;
for i=1:numFlies
    numTrials = length(allFlyData{i}.(trial));
    scatter(zeros(numTrials,1)+i,CCs(i,1:numTrials),'k');
    line([i-0.25 i+0.25], [mean(CCs(i,1:numTrials)) mean(CCs(i,1:numTrials))],'LineWidth',2,'Color','k');
end
xlim([0 numFlies+1]);

%Plot an example
flyID = 10;
darkID = 3;

% Get the unwrapped rotational position
Rot = UnWrap(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,2),2,0);
sgolayOrder = 3;
sgolayWindow = 11;
Rot = sgolayfilt(Rot,sgolayOrder,sgolayWindow);

% Find the unwrapped PVA
RSig = allFlyData{flyID}.(trial){darkID}.RROIaveMax;
GSig = allFlyData{flyID}.(trial){darkID}.GROIaveMax;

% Calculate the PVA for 60D05
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
angsraw = angsraw';
clear meanAngRaw;
clear meanIntRaw;
for ts = 1:length(RSig)
    meanAngRaw(ts) = circ_mean(angsraw, squeeze(RSig(:,ts)));
    meanIntRaw(ts) = circ_r(angsraw, squeeze(RSig(:,ts)));
end

PVA = UnWrap(meanAngRaw,1.5,0);

CC = corrcoef(PVA,Rot);
CCs(flyID,darkID) = CC(1,2);

exPVAPlot =  figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,1,1);
hold on;
imagesc([allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(1,1),allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(end,1)],...
    [-pi pi],...
    GSig-1);
axis tight;
set(gca,'FontSize',16);
caxis([0 2]);

subplot(4,1,2);
hold on;
imagesc([allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(1,1),allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(end,1)],...
    [-pi pi],...
    RSig-1);
colormap(brewermap(64, 'Blues'));
plot(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,1),meanAngRaw);
axis tight;
set(gca,'FontSize',16);
caxis([0 1]);

subplot(4,1,3);
hold on;
plot(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,1),Rot-Rot(1));
plot(allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(:,1),PVA-PVA(1));
axis tight;
set(gca,'FontSize',16);
xlabel('time (sec)');
ylabel('cumulatative rotation (rad)');
ylim([-2*pi 2*pi]);

subplot(4,1,4);
imagesc([allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(1,1),allFlyData{flyID}.(trial){darkID}.positionDatMatch.OffsetRotMatch(end,1)],...
    [-pi pi],...
    meanIntRaw);
caxis([0 0.2]);
colormap(brewermap(64, 'Greys'));
colorbar;

set(exPVAPlot,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(exPVAPlot,'PVAExampleWNOGreys','-dpdf');

%% Sort data into periods where the fly is walking forward with some velocity, but not turning

% Threshold velocities to consider
velThreshLCvF = [0:0.1:1.9];
velThreshUCvF = [0.1:0.1:2];
vRCut = 10000;
% velThreshLC = [0:45:135];
% velThreshUC = [45:45:180];
% XCVals = {};
% for flyID = 1:length(allFlyData)
%     XCVals{flyID}.ID = allFlyData{flyID}.ID;
% end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID})
            for velCut = 1:length(velThreshLCvF)
                % Sort into right and left turn bins by increasing or
                % decreasing velocity
                RBin = [];
                GBin = [];
                
                minFG = allFlyData{flyID}.dark{darkID}.positionDatMatch.minFG;
                maxFG = allFlyData{flyID}.dark{darkID}.positionDatMatch.maxFG;
                GROIaveMax = allFlyData{flyID}.dark{darkID}.GROIaveMax;
                RROIaveMax = allFlyData{flyID}.dark{darkID}.RROIaveMax;
                vFor = allFlyData{flyID}.dark{darkID}.positionDatMatch.vF;
                vR = allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot*180/pi;
                aFly = diff(vR);

                for i=1:length(GROIaveMax)-2
                    if vR(i) < vRCut 
                        if vFor(i) > velThreshLCvF(velCut) ...
                                & vFor(i) < velThreshUCvF(velCut)
                            RBin = horzcat(RBin, RROIaveMax(:,i));
                            GBin = horzcat(GBin, GROIaveMax(:,i));
                        end
                    end
                end

                % Save the values to the array
                XCVals{flyID}.velF{velCut}.threshold = velThreshLCvF(velCut);
                XCVals{flyID}.velF{velCut}.dark{darkID}.RBin = RBin;
                XCVals{flyID}.velF{velCut}.dark{darkID}.GBin = GBin;
            end
            
            clear RROIaveMax GROIaveMax vR vFor;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLCvF)       
       RBinGAllAlign = [];
       GBinGAllAlign = [];
       
       RBinRAllAlign = [];
       GBinRAllAlign = [];

        for trial=1:length(XCVals{flyID}.velF{velCut}.dark)
           RBinGAlign = [];
           GBinGAlign = [];

           RBinRAlign = [];
           GBinRAlign = [];
       
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.velF{velCut}.dark{trial})
                for i=1:size(XCVals{flyID}.velF{velCut}.dark{trial}.RBin,2)
                    pkPosRalign = find(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i)==...
                        max(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i)));
                    pkPosGalign = find(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i)==...
                        max(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i)));
                    RBinGAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i),num_ROIs/2-pkPosGalign);
                    GBinGAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i),num_ROIs/2-pkPosGalign);
                    RBinRAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i),num_ROIs/2-pkPosRalign);
                    GBinRAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i),num_ROIs/2-pkPosRalign);
                end
            end
            
            RBinGAllAlign = horzcat(RBinGAllAlign,RBinGAlign);
            GBinGAllAlign = horzcat(GBinGAllAlign,GBinGAlign);
            
            RBinRAllAlign = horzcat(RBinRAllAlign,RBinGAlign);
            GBinRAllAlign = horzcat(GBinRAllAlign,GBinGAlign);
            

       end
       XCVals{flyID}.velF{velCut}.RBinGAllAlign = RBinGAllAlign;
       XCVals{flyID}.velF{velCut}.GBinGAllAlign = GBinGAllAlign;
       XCVals{flyID}.velF{velCut}.RBinRAllAlign = RBinRAllAlign;
       XCVals{flyID}.velF{velCut}.GBinRAllAlign = GBinRAllAlign;
    end 
end

% figure;
% for flyID = 1:numFlies
%     for velCut = 1:length(velThreshLC)
%         subplot(numFlies,length(velThreshLC),velCut + length(velThreshLC)*(flyID-1))
%         hold on;
%         plot(mean(XCVals{flyID}.velF{velCut}.RBinRAllAlign,2),'r');
%         plot(mean(XCVals{flyID}.velF{velCut}.GBinRAllAlign,2),'g');
%     end
% end

RmaxVals = zeros(numFlies,length(velThreshLCvF));
GmaxVals = zeros(numFlies,length(velThreshLCvF));
for flyID = 1:numFlies
    for velCut = 1:length(velThreshLCvF)
        if ~isempty(XCVals{flyID}.velF{velCut}.RBinRAllAlign)
            
            RmaxVals(flyID,velCut) = max(mean(XCVals{flyID}.velF{velCut}.RBinGAllAlign,2))-min(mean(XCVals{flyID}.velF{velCut}.RBinGAllAlign,2));
            GmaxVals(flyID,velCut) = max(mean(XCVals{flyID}.velF{velCut}.GBinGAllAlign,2))-min(mean(XCVals{flyID}.velF{velCut}.GBinGAllAlign,2));
        end
    end
end

for flyID = 1:numFlies
    RnormVal = RmaxVals(flyID,1);
    GnormVal = GmaxVals(flyID,1);
    RmaxVals(flyID,:) = RmaxVals(flyID,:)./RnormVal;
    GmaxVals(flyID,:) = GmaxVals(flyID,:)./GnormVal;
end

Rmean = zeros(length(velCut),1);
Gmean = zeros(length(velCut),1);
for velCut = 1:length(velThreshLCvF)
    RNow = RmaxVals(:,velCut); 
    Rmean(velCut) = mean(RNow(find(RNow > 0)));
    GNow = GmaxVals(:,velCut); 
    Gmean(velCut) = mean(GNow(find(GNow > 0)));
end

figure;

for velCut = 1:length(velThreshLCvF)
    subplot(2,2,1)
    hold on;
    GPlot = GmaxVals(:,velCut);
    scatter(zeros(length(GPlot(find(GPlot~=0))),1)+velCut,GPlot(find(GPlot~=0)),'k','filled');
    line([velCut-0.5 velCut+0.5], [Gmean(velCut) Gmean(velCut)],'Color','k','LineWidth',2);
    set(gca,'FontSize',14);
    ylabel('Peak max');
    ylim([0.5 8]);
    xlim([0.5 (length(velThreshLCvF)+ 0.5)]);
    xlabel('vFor bins (cm/sec)');

    subplot(2,2,3)
    hold on;
    RPlot = RmaxVals(:,velCut);
    scatter(zeros(length(RPlot(find(RPlot~=0))),1)+velCut,RPlot(find(RPlot~=0)),'k','filled');
    line([velCut-0.5 velCut+0.5], [Rmean(velCut) Rmean(velCut)],'Color','k','LineWidth',2);
    ylim([0 4]);
    xlim([0.5 (length(velThreshLCvF) +0.5)]);
    set(gca,'FontSize',14);
    ylabel('Peak max');
    xlabel('vFor bins (cm/sec)');

end


%% Get the peak heights and widths according to rotational velocity
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=1:length(allFlyData)
    GallDat = [];
    RallDat = [];
    for trial = 1:length(allFlyData{flyID}.dark)
        GallDat = horzcat(GallDat,...
            reshape(allFlyData{flyID}.dark{trial}.GROIaveMax,...
            [1, size(allFlyData{flyID}.dark{trial}.GROIaveMax,1)*...
            size(allFlyData{flyID}.dark{trial}.GROIaveMax,2)]));
        RallDat = horzcat(RallDat,...
            reshape(allFlyData{flyID}.dark{trial}.RROIaveMax,...
            [1, size(allFlyData{flyID}.dark{trial}.RROIaveMax,1)*...
            size(allFlyData{flyID}.dark{trial}.RROIaveMax,2)]));
    end
    GSDNow = std(GallDat);
    RSDNow = std(RallDat);
    for velCut = 1:length(velThreshLC)     
       GturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec));
       GturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec));
       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec));
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec));

       if ~isempty(RturnsG) 
           % Calculate the PVA differences
           GpkMax(flyID, velCut) = mean(max(RturnsG,[],1)-min(RturnsG,[],1));
           GpkWids = zeros(size(RturnsG,2),1);
           
           for git = 1:size(RturnsG,2)
               XCNow = RturnsG(:,git);
               if max(XCNow) > mean(XCNow) + 0*GSDNow
                   overThresh = find((XCNow - min(XCNow))> 0.5*(max(XCNow)-min(XCNow)));
                   GpkWids(git) = length(overThresh)-1;
               end
           end
           GpkWid(flyID, velCut) = mean(GpkWids(find(GpkWids>0)));
               
           RpkMax(flyID, velCut) = mean(max(GturnsR,[],1)-min(GturnsR,[],1));
           RpkWids = zeros(size(GturnsR,2),1);
           for git = 1:size(GturnsR,2)
               XCNow = GturnsR(:,git);
               if max(XCNow) > mean(XCNow) + 0*RSDNow
                   overThresh = find((XCNow - min(XCNow))> 0.5*(max(XCNow)-min(XCNow)));
                   RpkWids(git) = length(overThresh)-1;
               end
           end
           RpkWid(flyID, velCut) = mean(RpkWids(find(RpkWids>0)));
       end
    end

end


GNorm = GpkMax(:, 1);
RNorm = RpkMax(:, 1);
for flyID = 1:size(GpkMax,1)
    for velCut = 1:length(velThreshLC)
        GpkMax(flyID,velCut) = GpkMax(flyID,velCut)/GNorm(flyID);
        RpkMax(flyID,velCut) = RpkMax(flyID,velCut)/RNorm(flyID);
    end
end
        

PVAdiff = PVAdiff;

figure;
set(gcf,'Position',[100 100 560*2.5 420*2]);
subplot(2,3,1);
hold on;
for i=1:length(velThreshLC)
   PltNow = GpkMax(:,i);
   scatter(zeros(length(PltNow(find(PltNow~=0))),1)+30*i-15,PltNow(find(PltNow~=0)),'filled','k'); 
   line([30*i-30 30*i],[mean(PltNow(find(PltNow~=0))) mean(PltNow(find(PltNow~=0)))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 180]);
ylim([0.75 2.5]);
xlabel('vr bins (deg/sec)');

subplot(2,3,3);
hold on;
for i=1:length(velThreshLC)
   PltNow = GpkWid(:,i);
   scatter(zeros(length(PltNow(find(PltNow~=0))),1)+30*i-15,22.5*PltNow(find(PltNow~=0)),'filled','k'); 
   line([30*i-15 30*i],[22.5*mean(PltNow(find(PltNow~=0))) 22.5*mean(PltNow(find(PltNow~=0)))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 180]);
xlabel('vr bins (deg/sec)');
ylim([0 150]);

subplot(2,3,4);
hold on;
for i=1:length(velThreshLC)
   PltNow = RpkMax(:,i);
   scatter(zeros(length(PltNow(find(PltNow~=0))),1)+30*i-15,PltNow(find(PltNow~=0)),'filled','k'); 
   line([30*i-30 30*i],[mean(PltNow(find(PltNow~=0))) mean(PltNow(find(PltNow~=0)))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 180]);
ylim([0.75 2]);
xlabel('vr bins (deg/sec)');

subplot(2,3,6);
hold on;
for i=1:length(velThreshLC)
   PltNow = RpkWid(:,i);
   scatter(zeros(length(PltNow(find(PltNow~=0))),1)+30*i-15,22.5*PltNow(find(PltNow~=0)),'filled','k'); 
   line([30*i-30 30*i],[22.5*mean(PltNow(find(PltNow~=0))) 22.5*mean(PltNow(find(PltNow~=0)))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 180]);
xlabel('vr bins (deg/sec)');
ylim([0 180]);

%% Get the peak heights and widths according to forward velocity

% Threshold velocities to consider
velThreshLCvF = [0:0.1:1.9];
velThreshUCvF = [0.1:0.1:2];
vRCut = 10000;

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID})
            for velCut = 1:length(velThreshLCvF)
                % Sort into right and left turn bins by increasing or
                % decreasing velocity
                RBin = [];
                GBin = [];
                
                minFG = allFlyData{flyID}.dark{darkID}.positionDatMatch.minFG;
                maxFG = allFlyData{flyID}.dark{darkID}.positionDatMatch.maxFG;
                GROIaveMax = allFlyData{flyID}.dark{darkID}.GROIaveMax;
                RROIaveMax = allFlyData{flyID}.dark{darkID}.RROIaveMax;
                vFor = allFlyData{flyID}.dark{darkID}.positionDatMatch.vF;
                vR = allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot*180/pi;
                aFly = diff(vR);

                for i=1:length(GROIaveMax)-2
                    if vR(i) < vRCut 
                        if vFor(i) > velThreshLCvF(velCut) ...
                                & vFor(i) < velThreshUCvF(velCut)
                            RBin = horzcat(RBin, RROIaveMax(:,i));
                            GBin = horzcat(GBin, GROIaveMax(:,i));
                        end
                    end
                end

                % Save the values to the array
                XCVals{flyID}.velF{velCut}.threshold = velThreshLCvF(velCut);
                XCVals{flyID}.velF{velCut}.dark{darkID}.RBin = RBin;
                XCVals{flyID}.velF{velCut}.dark{darkID}.GBin = GBin;
            end
            
            clear RROIaveMax GROIaveMax vR vFor;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLCvF)       
       RBinGAllAlign = [];
       GBinGAllAlign = [];
       
       RBinRAllAlign = [];
       GBinRAllAlign = [];

        for trial=1:length(XCVals{flyID}.velF{velCut}.dark)
           RBinGAlign = [];
           GBinGAlign = [];

           RBinRAlign = [];
           GBinRAlign = [];
       
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.velF{velCut}.dark{trial})
                for i=1:size(XCVals{flyID}.velF{velCut}.dark{trial}.RBin,2)
                    pkPosRalign = find(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i)==...
                        max(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i)));
                    pkPosGalign = find(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i)==...
                        max(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i)));
                    RBinGAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i),num_ROIs/2-pkPosGalign);
                    GBinGAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i),num_ROIs/2-pkPosGalign);
                    RBinRAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.RBin(:,i),num_ROIs/2-pkPosRalign);
                    GBinRAlign(:,i) = circshift(XCVals{flyID}.velF{velCut}.dark{trial}.GBin(:,i),num_ROIs/2-pkPosRalign);
                end
            end
            
            RBinGAllAlign = horzcat(RBinGAllAlign,RBinGAlign);
            GBinGAllAlign = horzcat(GBinGAllAlign,GBinGAlign);
            
            RBinRAllAlign = horzcat(RBinRAllAlign,RBinGAlign);
            GBinRAllAlign = horzcat(GBinRAllAlign,GBinGAlign);
            

       end
       XCVals{flyID}.velF{velCut}.RBinGAllAlign = RBinGAllAlign;
       XCVals{flyID}.velF{velCut}.GBinGAllAlign = GBinGAllAlign;
       XCVals{flyID}.velF{velCut}.RBinRAllAlign = RBinRAllAlign;
       XCVals{flyID}.velF{velCut}.GBinRAllAlign = GBinRAllAlign;
    end 
end

% figure;
% for flyID = 1:numFlies
%     for velCut = 1:length(velThreshLC)
%         subplot(numFlies,length(velThreshLC),velCut + length(velThreshLC)*(flyID-1))
%         hold on;
%         plot(mean(XCVals{flyID}.velF{velCut}.RBinRAllAlign,2),'r');
%         plot(mean(XCVals{flyID}.velF{velCut}.GBinRAllAlign,2),'g');
%     end
% end

RmaxVals = zeros(numFlies,length(velThreshLCvF));
GmaxVals = zeros(numFlies,length(velThreshLCvF));
for flyID = 1:numFlies
    for velCut = 1:length(velThreshLCvF)
        if ~isempty(XCVals{flyID}.velF{velCut}.RBinRAllAlign)
            
            RmaxVals(flyID,velCut) = mean(max(XCVals{flyID}.velF{velCut}.RBinGAllAlign,[],1)-min(XCVals{flyID}.velF{velCut}.RBinGAllAlign,[],1));
            GmaxVals(flyID,velCut) = mean(max(XCVals{flyID}.velF{velCut}.GBinGAllAlign,[],1)-min(XCVals{flyID}.velF{velCut}.GBinGAllAlign,[],1));
        end
    end
end

for flyID = 1:numFlies
    RnormVal = RmaxVals(flyID,1);
    GnormVal = GmaxVals(flyID,1);
    RmaxVals(flyID,:) = RmaxVals(flyID,:)./RnormVal;
    GmaxVals(flyID,:) = GmaxVals(flyID,:)./GnormVal;
end

Rmean = zeros(length(velCut),1);
Gmean = zeros(length(velCut),1);
for velCut = 1:length(velThreshLCvF)
    RNow = RmaxVals(:,velCut); 
    Rmean(velCut) = mean(RNow(find(RNow > 0)));
    GNow = GmaxVals(:,velCut); 
    Gmean(velCut) = mean(GNow(find(GNow > 0)));
end

for velCut = 1:length(velThreshLCvF)
    subplot(2,3,2)
    hold on;
    GPlot = GmaxVals(:,velCut);
    scatter(zeros(length(GPlot(find(GPlot~=0))),1)+velCut*0.1-0.05,GPlot(find(GPlot~=0)),'k','filled');
    line([0.1*velCut-0.1 0.1*velCut], [Gmean(velCut) Gmean(velCut)],'Color','k','LineWidth',2);
    set(gca,'FontSize',14);
    ylabel('Peak max');
    ylim([0.75 2.5]);
    xlim([0 0.6]);
    xlabel('vFor bins (cm/sec)');

    subplot(2,3,5)
    hold on;
    RPlot = RmaxVals(:,velCut);
    scatter(zeros(length(RPlot(find(RPlot~=0))),1)+velCut*0.1-0.05,RPlot(find(RPlot~=0)),'k','filled');
    line([0.1*velCut-0.1 0.1*velCut], [Rmean(velCut) Rmean(velCut)],'Color','k','LineWidth',2);
    ylim([0.75 2]);
    xlim([0 0.6]);
    set(gca,'FontSize',14);
    ylabel('Peak max');
    xlabel('vFor bins (cm/sec)');

end
%% Load all of the data
% Specify directories
dirRoot = 'D:\Imaging\2Color';

imDir{1} = '20160518';
imDir{2} = '20160519';
redLine = '60D05';
greenLine = '37F06';


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
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
            for flyID = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyID}.ID)
                    allFlyData{flyID}.openLoop{str2num(fileName(end-4))}.positionDat = positionDat;
                    allFlyData{flyID}.openLoop{str2num(fileName(end-4))}.RROIaveMax = RROIaveMax;
                    allFlyData{flyID}.openLoop{str2num(fileName(end-4))}.GROIaveMax = GROIaveMax;
                    allFlyData{flyID}.openLoop{str2num(fileName(end-4))}.RROIaveMean = RROIaveMean;
                    allFlyData{flyID}.openLoop{str2num(fileName(end-4))}.GROIaveMean = GROIaveMean;
                end
            end
        end 
    end
end

%% Compare the bump shifts for different velocity thresholds for the openLoop trials
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
velThreshLC = [0 20 40 60];
velThreshUC = [20 40 60 80];
XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the openLoop cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the openLoop trials
    for openLoopID = 1:length(allFlyData{flyID}.openLoop)
        if ~isempty(allFlyData{flyID}.openLoop{openLoopID})
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.openLoop{openLoopID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.openLoop{openLoopID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.openLoop{openLoopID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.openLoop{openLoopID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.openLoop{openLoopID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.openLoop{openLoopID}.RROIaveMax));
            tSpan = allFlyData{flyID}.openLoop{openLoopID}.positionDat.t - ...
                allFlyData{flyID}.openLoop{openLoopID}.positionDat.t(1) + ...
                allFlyData{flyID}.openLoop{openLoopID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.openLoop{openLoopID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.openLoop{openLoopID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.openLoop{openLoopID}.positionDat.tFrameGrab)/num_planes);

            % Match the position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            OffsetTransMatch = zeros(maxFG-minFG+1,1);
            OffsetDirMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.openLoop{openLoopID}.positionDat.t >=...
                    (allFlyData{flyID}.openLoop{openLoopID}.positionDat.t(1) +...
                    (allFlyData{flyID}.openLoop{openLoopID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.openLoop{openLoopID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.openLoop{openLoopID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.openLoop{openLoopID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.openLoop{openLoopID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.openLoop{openLoopID}.positionDat.OffsetLat(tMatch(1));
                OffsetTransMatch(interp-minFG+1) = allFlyData{flyID}.openLoop{openLoopID}.positionDat.trans(tMatch(1));
                OffsetDirMatch(interp-minFG+1) = allFlyData{flyID}.openLoop{openLoopID}.positionDat.direction(tMatch(1));
            end
            
            % Unwrap the rotation and smooth over 1 sec to get periods of
            % continuous turning            
            OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),1.5,0);
            
            OffsetRotCorrected = OffsetRotUnwrap;
            flip = 0;
            for crt = 1:length(OffsetRotCorrected)-1
                if flip==0
                    if OffsetTransMatch(crt+1) == 0 & OffsetTransMatch(crt) == 2
                        OffsetRotCorrected(crt+1:end) = OffsetRotCorrected(crt+1:end) - (OffsetRotMatch(crt+1,2)-OffsetRotMatch(crt,2));
                        flip = 1;
                    end
                else
                    if OffsetTransMatch(crt+1) == 2 & OffsetTransMatch(crt) == 0
                        OffsetRotCorrected(crt+1:end) = OffsetRotCorrected(crt+1:end) - (OffsetRotMatch(crt+1,2)-OffsetRotMatch(crt,2));
                        flip = 0;
                    end
                end
            end          
                        
            OffsetRotUnwrap = smooth(OffsetRotUnwrap,round(tSmooth./mean(diff(OffsetRotMatch(:,1)))));
            stripeInView = mod(OffsetRotUnwrap+pi,2*pi)-pi;
            
            [posRot, posFor, posLat] = PositionConverter(allFlyData{flyID}.openLoop{openLoopID}.positionDat);
            
            % Match the fly position data to the framegrab times
            posRotMatch = zeros(maxFG-minFG+1,1);
            posForMatch = zeros(maxFG-minFG+1,1);
            posLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.openLoop{openLoopID}.positionDat.t >=...
                    (allFlyData{flyID}.openLoop{openLoopID}.positionDat.t(1) +...
                    (allFlyData{flyID}.openLoop{openLoopID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-positionDat.tVR(1))/10000));
                posRotMatch(interp-minFG+1) = pi/180*posRot(tMatch(1));
                posForMatch(interp-minFG+1) = posFor(tMatch(1));
                posLatMatch(interp-minFG+1) = posLat(tMatch(1));
            end
            
            posRotUnwrap = UnWrap(posRotMatch,1.5,0);
            
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
                    if posRotUnwrap(i+1-minFG) - posRotUnwrap(i-minFG) > velThreshLCNow ...
                            & posRotUnwrap(i+1-minFG) - posRotUnwrap(i-minFG) < velThreshUCNow ...
                            & OffsetTransMatch(i+1-minFG) == 2 & abs(stripeInView(i+1-minFG)) < 2*pi/3 & OffsetDirMatch(i+1-minFG) == -1
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        RturnsG = horzcat(RturnsG, GROIaveMaxFilt(:,i));
                        RturnsR = horzcat(RturnsR, RROIaveMaxFilt(:,i));
                    elseif posRotUnwrap(i+1-minFG) - posRotUnwrap(i-minFG) < -velThreshLCNow ...
                            & posRotUnwrap(i+1-minFG) - posRotUnwrap(i-minFG) > -velThreshUCNow ...
                            & OffsetTransMatch(i+1-minFG) == 2 & abs(stripeInView(i+1-minFG)) < 2*pi/3 & OffsetDirMatch(i+1-minFG) == -1
%                             & meanIntRawG(i) > 0.5*max(meanIntRawG) & meanIntRawR(i) > 0.5*max(meanIntRawR)
                        LturnsG = horzcat(LturnsG, GROIaveMaxFilt(:,i));
                        LturnsR = horzcat(LturnsR, RROIaveMaxFilt(:,i));
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.openLoop{openLoopID}.RturnsG = RturnsG;
                XCVals{flyID}.vel{velCut}.openLoop{openLoopID}.RturnsR = RturnsR;
                XCVals{flyID}.vel{velCut}.openLoop{openLoopID}.LturnsG = LturnsG;
                XCVals{flyID}.vel{velCut}.openLoop{openLoopID}.LturnsR = LturnsR;
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
       for trial=1:length(XCVals{flyID}.vel{velCut}.openLoop)
            RturnsGRalign = [];
            RturnsGGalign = [];
            RturnsRRalign = [];
            RturnsRGalign = [];
            LturnsGRalign = [];
            LturnsGGalign = [];
            LturnsRRalign = [];
            LturnsRGalign = [];
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.openLoop{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsG(:,i)));
                    RturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsG(:,i),num_ROIs/2-pkPosRalign);
                    RturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsG(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsR(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.RturnsR(:,i),num_ROIs/2-pkPosGalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsR,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsR(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsR(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsG(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsG(:,i)));
                    LturnsGRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsG(:,i),num_ROIs/2-pkPosRalign);
                    LturnsGGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsG(:,i),num_ROIs/2-pkPosGalign);
                    LturnsRRalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsR(:,i),num_ROIs/2-pkPosRalign);
                    LturnsRGalign(:,i) = circshift(XCVals{flyID}.vel{velCut}.openLoop{trial}.LturnsR(:,i),num_ROIs/2-pkPosGalign);
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
       XCVals{flyID}.vel{velCut}.RturnsGAllopenLoopRalign = RturnsGallRalign;
       XCVals{flyID}.vel{velCut}.RturnsGAllopenLoopGalign = RturnsGallGalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllopenLoopRalign = RturnsRallRalign;
       XCVals{flyID}.vel{velCut}.RturnsRAllopenLoopGalign = RturnsRallGalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllopenLoopRalign = LturnsGallRalign;
       XCVals{flyID}.vel{velCut}.LturnsGAllopenLoopGalign = LturnsGallGalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllopenLoopRalign = LturnsRallRalign;
       XCVals{flyID}.vel{velCut}.LturnsRAllopenLoopGalign = LturnsRallGalign;
    end 
end

bumpCompG = figure('units','normalized','outerposition',[0 0 0.5 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.openLoop)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllopenLoopRalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllopenLoopRalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllopenLoopRalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllopenLoopRalign;
           if ~isempty(RturnsG) 
               RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
               RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           else
               RnetR = [];
               RnetG = [];
           end
           if ~isempty(LturnsG)
               LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
               LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           else
               LnetR = [];
               LnetG = [];
           end
           
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
       for trial=1:length(XCVals{flyID}.vel{velCut}.openLoop)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllopenLoopGalign;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllopenLoopGalign;
           LturnsG = XCVals{flyID}.vel{velCut}.LturnsGAllopenLoopGalign;
           LturnsR = XCVals{flyID}.vel{velCut}.LturnsRAllopenLoopGalign;
           if ~isempty(RturnsG) 
               RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
               RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           else
               RnetR = [];
               RnetG = [];
           end
           if ~isempty(LturnsG)
               LnetR = (LturnsR-min(median(LturnsR,2)))./(max(median(LturnsR,2))-min(median(LturnsR,2)));
               LnetG = (LturnsG-min(median(LturnsG,2)))./(max(median(LturnsG,2))-min(median(LturnsG,2)));
           else
               LnetR = [];
               LnetG = [];
           end
           
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
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
LPVAdiffMean = zeros(numFlies,plotRange*2+1);
LPeakdiff = zeros(numFlies,plotRange*2+1);
LPeakdiffMean = zeros(numFlies,plotRange*2+1);
RPVAdiff = zeros(numFlies,plotRange*2+1);
RPVAdiffMean = zeros(numFlies,plotRange*2+1);
RPeakdiff = zeros(numFlies,plotRange*2+1);
RPeakdiffMean = zeros(numFlies,plotRange*2+1);

for flyID = 1:length(allFlyData)
    
    LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
    RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
    LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
    RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;
    
    % Calculate the PVA difference snd peak difference of the mean
    LBinG = mean(LPBDataGAllStop,2);
    RBinG = mean(RPBDataGAllStop,2);
    LBinR = mean(LPBDataRAllStop,2);
    RBinR = mean(RPBDataRAllStop,2);
    
    LGMax = max(LBinG);
    RGMax = max(RBinG);
    LRMax = max(LBinR);
    RRMax = max(RBinR);
    
    LPVAdiff(flyID,plotRange+1) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
    LPVAdiff(flyID,plotRange+1) = min(LPVAdiff(flyID,plotRange+1),9-LPVAdiff(flyID,plotRange+1));
    LPeakdiff(flyID,plotRange+1) = greenSpan(1) + find(LBinG==LGMax)-...
        (redSpan(1) + find(LBinR ==LRMax));
    LPeakdiff(flyID,plotRange+1) = min(LPeakdiff(flyID,plotRange+1),9-LPeakdiff(flyID,plotRange+1));

    RPVAdiff(flyID,plotRange+1) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
    RPVAdiff(flyID,plotRange+1) = min(RPVAdiff(flyID,plotRange+1),9-RPVAdiff(flyID,plotRange+1));
    RPeakdiff(flyID,plotRange+1) = greenSpan(1) + find(RBinG==RGMax)-...
        (redSpan(1) + find(RBinR ==RRMax));
    RPeakdiff(flyID,plotRange+1) = min(RPeakdiff(flyID,plotRange+1),9-RPeakdiff(flyID,plotRange+1));
    
    % Calculate the mean of the individual PVA and peak differences
    LPVAdiffMeanTemp = zeros(1,size(LPBDataRAllStop,2));
    LPeakdiffMeanTemp = zeros(1,size(LPBDataRAllStop,2));
    RPVAdiffMeanTemp = zeros(1,size(RPBDataRAllStop,2));
    RPeakdiffMeanTemp = zeros(1,size(RPBDataRAllStop,2));
    
    for i=1:size(LPBDataRAllStop,2)
       LPVAdiffMeanTemp(i) =  greenSpan(1) + circ_mean(angsraw,LPBDataGAllStop(greenSpan(1:8),i))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LPBDataRAllStop(redSpan(1:8),i))*num_ROIs/(2*pi));
       RPVAdiffMeanTemp(i) = -greenSpan(1) + circ_mean(angsraw,RPBDataGAllStop(sort(10-greenSpan(1:8)),i))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RPBDataRAllStop(sort(10-redSpan(1:8)),i))*num_ROIs/(2*pi));
    
       LGMaxTemp = max(LPBDataGAllStop(greenSpan(1:8),i));
       RGMaxTemp = max(RPBDataGAllStop(sort(10-greenSpan(1:8)),i));
       LRMaxTemp = max(LPBDataRAllStop(redSpan(1:8),i));
       RRMaxTemp = max(RPBDataRAllStop(sort(10-redSpan(1:8)),i));
       
       LPeakdiffMeanTemp(i) =  greenSpan(1) + find(LPBDataGAllStop(greenSpan(1:8),i)==LGMaxTemp) -...
        (redSpan(1) + find(LPBDataRAllStop(redSpan(1:8),i)==LRMaxTemp));
       RPeakdiffMeanTemp(i) = -greenSpan(1) + find(RPBDataGAllStop(sort(10-greenSpan(1:8)),i)==RGMaxTemp)-...
        (-redSpan(1) + find(RPBDataRAllStop(sort(10-redSpan(1:8)),i)==RRMaxTemp));
    end
    LPVAdiffMean(flyID,plotRange+1) = mean(LPVAdiffMeanTemp);
    RPVAdiffMean(flyID,plotRange+1) = mean(RPVAdiffMeanTemp);
    LPeakdiffMean(flyID,plotRange+1) = mean(LPeakdiffMeanTemp);
    RPeakdiffMean(flyID,plotRange+1) = mean(RPeakdiffMeanTemp);
    
    
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
            LGMax = max(LBinG);
            RGMax = max(RBinG);
            LRMax = max(LBinR);
            RRMax = max(RBinR);
    
            LPVAdiff(flyID,plotRange+1-binID) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
            %LPVAdiff(flyID,plotRange+1-binID) = min(LPVAdiff(flyID,plotRange+1-binID),9-LPVAdiff(flyID,plotRange+1-binID));
            LPeakdiff(flyID,plotRange+1-binID) = greenSpan(1) + find(LBinG==LGMax)-...
            (redSpan(1) + find(LBinR ==LRMax));
%             LPeakdiff(flyID,plotRange+1-binID) = min(LPeakdiff(flyID,plotRange+1),9-LPeakdiff(flyID,plotRange+1));

            RPVAdiff(flyID,plotRange+1-binID) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
            %RPVAdiff(flyID,plotRange+1-binID) = min(RPVAdiff(flyID,plotRange+1-binID),9-RPVAdiff(flyID,plotRange+1-binID));
            RPeakdiff(flyID,plotRange+1-binID) = greenSpan(1) + find(RBinG==RGMax)-...
            (redSpan(1) + find(RBinR ==RRMax));
%             RPeakdiff(flyID,plotRange+1-binID) = min(RPeakdiff(flyID,plotRange+1),9-RPeakdiff(flyID,plotRange+1));

            % Calculate the mean of the individual PVA and peak differences
            LPVAdiffMeanTemp = zeros(1,size(LPBDataRAllCW,2));
            LPeakdiffMeanTemp = zeros(1,size(LPBDataRAllCW,2));
            RPVAdiffMeanTemp = zeros(1,size(RPBDataRAllCW,2));
            RPeakdiffMeanTemp = zeros(1,size(RPBDataRAllCW,2));

            for i=1:size(LPBDataRAllCW,2)
               LPVAdiffMeanTemp(i) =  greenSpan(1) + circ_mean(angsraw,LPBDataGAllCW(greenSpan(1:8),i))*num_ROIs/(2*pi)-...
                (redSpan(1) + circ_mean(angsraw,LPBDataRAllCW(redSpan(1:8),i))*num_ROIs/(2*pi));
               RPVAdiffMeanTemp(i) = -greenSpan(1) + circ_mean(angsraw,RPBDataGAllCW(sort(10-greenSpan(1:8)),i))*num_ROIs/(2*pi)-...
                (-redSpan(1) + circ_mean(angsraw,RPBDataRAllCW(sort(10-redSpan(1:8)),i))*num_ROIs/(2*pi));

               LGMaxTemp = max(LPBDataGAllCW(greenSpan(1:8),i));
               RGMaxTemp = max(RPBDataGAllCW(sort(10-greenSpan(1:8)),i));
               LRMaxTemp = max(LPBDataRAllCW(redSpan(1:8),i));
               RRMaxTemp = max(RPBDataRAllCW(sort(10-redSpan(1:8)),i));

               LPeakdiffMeanTemp(i) =  greenSpan(1) + find(LPBDataGAllCW(greenSpan(1:8),i)==LGMaxTemp) -...
                (redSpan(1) + find(LPBDataRAllCW(redSpan(1:8),i)==LRMaxTemp));
               RPeakdiffMeanTemp(i) = -greenSpan(1) + find(RPBDataGAllCW(sort(10-greenSpan(1:8)),i)==RGMaxTemp)-...
                (-redSpan(1) + find(RPBDataRAllCW(sort(10-redSpan(1:8)),i)==RRMaxTemp));
            end
            LPVAdiffMean(flyID,plotRange+1-binID) = mean(LPVAdiffMeanTemp);
            RPVAdiffMean(flyID,plotRange+1-binID) = mean(RPVAdiffMeanTemp);
            LPeakdiffMean(flyID,plotRange+1-binID) = mean(LPeakdiffMeanTemp);
            RPeakdiffMean(flyID,plotRange+1-binID) = mean(RPeakdiffMeanTemp);
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
            LGMax = max(LBinG);
            RGMax = max(RBinG);
            LRMax = max(LBinR);
            RRMax = max(RBinR);
            LPVAdiff(flyID,plotRange+binID) = greenSpan(1) + circ_mean(angsraw,LBinG(greenSpan(1:8)))*num_ROIs/(2*pi)-...
        (redSpan(1) + circ_mean(angsraw,LBinR(redSpan(1:8)))*num_ROIs/(2*pi));
%             LPVAdiff(flyID,plotRange+binID) = min(LPVAdiff(flyID,plotRange+binID),9-LPVAdiff(flyID,plotRange+binID));
            LPeakdiff(flyID,plotRange+1+binID) = greenSpan(1) + find(LBinG==LGMax)-...
            (redSpan(1) + find(LBinR ==LRMax));
%             LPeakdiff(flyID,plotRange+1+binID) = min(LPeakdiff(flyID,plotRange+1),9-LPeakdiff(flyID,plotRange+1));

            RPVAdiff(flyID,plotRange+binID) = -greenSpan(1) + circ_mean(angsraw,RBinG(sort(10-greenSpan(1:8))))*num_ROIs/(2*pi)-...
        (-redSpan(1) + circ_mean(angsraw,RBinR(sort(10-redSpan(1:8))))*num_ROIs/(2*pi));
%             RPVAdiff(flyID,plotRange+binID) = min(RPVAdiff(flyID,plotRange+binID),9-RPVAdiff(flyID,plotRange+binID));
            RPeakdiff(flyID,plotRange+1+binID) = greenSpan(1) + find(RBinG==RGMax)-...
            (redSpan(1) + find(RBinR ==RRMax));
%             RPeakdiff(flyID,plotRange+1+binID) = min(RPeakdiff(flyID,plotRange+1),9-RPeakdiff(flyID,plotRange+1));


            % Calculate the mean of the individual PVA and peak differences
            LPVAdiffMeanTemp = zeros(1,size(LPBDataRAllCCW,2));
            LPeakdiffMeanTemp = zeros(1,size(LPBDataRAllCCW,2));
            RPVAdiffMeanTemp = zeros(1,size(RPBDataRAllCCW,2));
            RPeakdiffMeanTemp = zeros(1,size(RPBDataRAllCCW,2));

            for i=1:size(LPBDataRAllCCW,2)
               LPVAdiffMeanTemp(i) =  greenSpan(1) + circ_mean(angsraw,LPBDataGAllCCW(greenSpan(1:8),i))*num_ROIs/(2*pi)-...
                (redSpan(1) + circ_mean(angsraw,LPBDataRAllCCW(redSpan(1:8),i))*num_ROIs/(2*pi));
               RPVAdiffMeanTemp(i) = -greenSpan(1) + circ_mean(angsraw,RPBDataGAllCCW(sort(10-greenSpan(1:8)),i))*num_ROIs/(2*pi)-...
                (-redSpan(1) + circ_mean(angsraw,RPBDataRAllCCW(sort(10-redSpan(1:8)),i))*num_ROIs/(2*pi));

               LGMaxTemp = max(LPBDataGAllCCW(greenSpan(1:8),i));
               RGMaxTemp = max(RPBDataGAllCCW(sort(10-greenSpan(1:8)),i));
               LRMaxTemp = max(LPBDataRAllCCW(redSpan(1:8),i));
               RRMaxTemp = max(RPBDataRAllCCW(sort(10-redSpan(1:8)),i));

               LPeakdiffMeanTemp(i) =  greenSpan(1) + find(LPBDataGAllCCW(greenSpan(1:8),i)==LGMaxTemp) -...
                (redSpan(1) + find(LPBDataRAllCCW(redSpan(1:8),i)==LRMaxTemp));
               RPeakdiffMeanTemp(i) = -greenSpan(1) + find(RPBDataGAllCCW(sort(10-greenSpan(1:8)),i)==RGMaxTemp)-...
                (-redSpan(1) + find(RPBDataRAllCCW(sort(10-redSpan(1:8)),i)==RRMaxTemp));
            end
            LPVAdiffMean(flyID,plotRange+1+binID) = mean(LPVAdiffMeanTemp);
            RPVAdiffMean(flyID,plotRange+1+binID) = mean(RPVAdiffMeanTemp);
            LPeakdiffMean(flyID,plotRange+1+binID) = mean(LPeakdiffMeanTemp);
            RPeakdiffMean(flyID,plotRange+1+binID) = mean(RPeakdiffMeanTemp);
        end
    end
    
end

PBFig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LPVAdiff(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RPVAdiff(:,velPlt),'r');
end
ylim([-2 2]);
subplot(2,2,2);
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LPeakdiff(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RPeakdiff(:,velPlt),'r');
end
ylim([-2 2]);
subplot(2,2,3);
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LPVAdiffMean(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RPVAdiffMean(:,velPlt),'r');
end
ylim([-2 2]);
subplot(2,2,4);
hold on;
for velPlt = 1:2*plotRange+1
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),LPeakdiffMean(:,velPlt),'g');
    scatter(zeros(numFlies,1)+vRSpan*(velPlt-1-plotRange),RPeakdiffMean(:,velPlt),'r');
end
ylim([-2 2]);

save(strcat('PVADiff_',greenLine,'-',redLine,'2.mat'),'LPVAdiff','RPVAdiff',...
    'LPeakdiff','RPeakdiff','LPVAdiffMean','RPVAdiffMean','LPeakdiffMean','RPeakdiffMean');

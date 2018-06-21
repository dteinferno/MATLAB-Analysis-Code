function [dataR, dataG] = BumpAlign(allFlyData, Min, Max, Span, RG, LR, trial)
    
%function to sort data into bins

%input:
%allFlyData = the data to be sorted
%Min = the minimum velocity bin (deg)
%Max = the maximum velocity bin (deg)
%Span = the span of each bin (deg)
%RG = Switch to tell whether 60D05 is green (0) or red (1)
%LR = Switch to tell whether we're looking at the left PB (0) or right
%PB (1)
%trial = the trial to be considered

%output:
%dataR = the sorted data for the red channel
%dataG = the sorted data for the green channel

% Specify the range of rotational velocities to consider
vRMin = Min*pi/180;
vRMax = Max*pi/180;
vRSpan = Span*pi/180;
vRBinNum = round((vRMax-vRMin)/vRSpan);

% Initilaize the sorted output data arrays
dataR = cell(size(allFlyData,2));
dataG = cell(size(allFlyData,2));

% Specify which range of ROIs to consider
if LR
    if RG
        rangeR = [10:17];
        rangeG = [11:18];
        rangeAll = [10:18];
    else
        rangeR = [11:18];
        rangeG = [10:17];
        rangeAll = [10:18];
    end
else
    if RG
        rangeR = [2:9];
        rangeG = [1:8];
        rangeAll = [1:9];
    else
        rangeR = [1:8];
        rangeG = [2:9];
        rangeAll = [1:9];
    end
end

% Sort the data
for flyID = 1:length(allFlyData)
    dataR{flyID}.CW = cell(vRBinNum,1);
    dataR{flyID}.CCW = cell(vRBinNum,1);
    dataR{flyID}.Stop = [];  
    
    dataG{flyID}.CW = cell(vRBinNum,1);
    dataG{flyID}.CCW = cell(vRBinNum,1);
    dataG{flyID}.Stop = [];  

    for trialID = 1:length(allFlyData{flyID}.(trial))
        vRot = allFlyData{flyID}.(trial){trialID}.positionDatMatch.vRot;
        RSig = allFlyData{flyID}.(trial){trialID}.RROIaveMax;
        GSig = allFlyData{flyID}.(trial){trialID}.GROIaveMax;
        for tStep = 1:length(vRot)
            if RG
                pkMax = max(GSig(rangeG,tStep+1));
                pkShift = -find(GSig(rangeG,tStep+1) == pkMax);
            else
                pkMax = max(RSig(rangeR,tStep+1));
                pkShift = -find(RSig(rangeR,tStep+1) == pkMax);
            end
            RSigNow = circshift(RSig(rangeAll,tStep+1),pkShift);
            GSigNow = circshift(GSig(rangeAll,tStep+1),pkShift);
            
            if vRot(tStep) == 0
                dataR{flyID}.Stop = horzcat(dataR{flyID}.Stop, RSigNow);
                dataG{flyID}.Stop = horzcat(dataG{flyID}.Stop, GSigNow);

            elseif vRot(tStep) > 0
                dataR{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataR{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    RSigNow);
                dataG{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataG{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    GSigNow);
            else
                dataR{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataR{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    RSigNow);
                dataG{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataG{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    GSigNow);
            end
        end
    end
end
    
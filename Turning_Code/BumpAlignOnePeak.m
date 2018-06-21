function [dataRedLPB, dataRedRPB, dataGreenLPB, dataGreenRPB] = BumpAlignOnePeak(allFlyData, Min, Max, Span, rangeR, rangeG, trial, glomShift, RGAlign)
    
%function to register data to one peak and sort into bins by roational velocity

%input:
%allFlyData = the data to be sorted
%Min = the minimum velocity bin (deg)
%Max = the maximum velocity bin (deg)
%Span = the span of each bin (deg)
%rangeR = glmeruli to consider for the red channel
%rangeG = glmeruli to consider for the green channel
%trial = the trial to be considered
%glomShift = number of glomeruli to shift the aligned peak
%RGAlign = 0 to align to the red peak and 1 to align to the green peak

%output:
%dataRedLPB = the sorted data for the left PB red channel
%dataRedRPB = the sorted data for the right PB red channel
%dataGreenLPB = the sorted data for the left PB green channel
%dataGreenRPB = the sorted data for the right PB green channel

% Specify the range of rotational velocities to consider
vRMin = Min*pi/180;
vRMax = Max*pi/180;
vRSpan = Span*pi/180;
vRBinNum = round((vRMax-vRMin)/vRSpan);

% Initilaize the sorted output data arrays
dataRedLPB = cell(size(allFlyData,2));
dataRedRPB = cell(size(allFlyData,2));
dataGreenLPB = cell(size(allFlyData,2));
dataGreenRPB = cell(size(allFlyData,2));

% Register and sort the data
for flyID = 1:length(allFlyData)
    
    % Initialize the cell arrays 
    dataRedLPB{flyID}.CW = cell(vRBinNum,1);
    dataRedLPB{flyID}.CCW = cell(vRBinNum,1);
    dataRedLPB{flyID}.Stop = [];  
    dataRedRPB{flyID}.CW = cell(vRBinNum,1);
    dataRedRPB{flyID}.CCW = cell(vRBinNum,1);
    dataRedRPB{flyID}.Stop = [];  
    
    dataGreenLPB{flyID}.CW = cell(vRBinNum,1);
    dataGreenLPB{flyID}.CCW = cell(vRBinNum,1);
    dataGreenLPB{flyID}.Stop = [];  
    dataGreenRPB{flyID}.CW = cell(vRBinNum,1);
    dataGreenRPB{flyID}.CCW = cell(vRBinNum,1);
    dataGreenRPB{flyID}.Stop = [];  

    % Step through the trials
    for trialID = 1:length(allFlyData{flyID}.(trial))
        
        % Load the data
        vRot = allFlyData{flyID}.(trial){trialID}.positionDatMatch.vRot;
        RSig = allFlyData{flyID}.(trial){trialID}.RROIaveMax;
        GSig = allFlyData{flyID}.(trial){trialID}.GROIaveMax;
        
        % Register and sort at each time point
        for tStep = 1:length(vRot)
            
            % Specify which color channel to use for peak alignment
            if RGAlign == 1
                pkMax = max(GSig(rangeG,tStep+1));
                pkShift = -find(GSig(rangeG,tStep+1) == pkMax);
            elseif RGAlign == 0
                pkMax = max(RSig(rangeR,tStep+1));
                pkShift = -find(RSig(rangeR,tStep+1) == pkMax);
            else
                disp('Error');
                break;
            end
            
            % Initialize arrays to hold the registered data
            RSigNowLPB = zeros(9,1);
            RSigNowRPB = zeros(9,1);
            GSigNowLPB = zeros(9,1);
            GSigNowRPB = zeros(9,1);
            
            % Register the data
            if rangeR(1) < 9
                RSigNowLPB(rangeR) = circshift(RSig(rangeR,tStep+1),pkShift+glomShift);
                RSigNowRPB(rangeG) = circshift(RSig(sort(19-rangeR),tStep+1),pkShift+glomShift);
                GSigNowLPB(rangeG) = circshift(GSig(rangeG,tStep+1),pkShift+glomShift);
                GSigNowRPB(rangeR) = circshift(GSig(sort(19-rangeG),tStep+1),pkShift+glomShift);
            else
                RSigNowLPB(sort(19-rangeR)) = circshift(RSig(sort(19-rangeR),tStep+1),pkShift+glomShift);
                RSigNowRPB(sort(19-rangeG)) = circshift(RSig(rangeR,tStep+1),pkShift+glomShift);
                GSigNowLPB(sort(19-rangeG)) = circshift(GSig(sort(19-rangeG),tStep+1),pkShift+glomShift);
                GSigNowRPB(sort(19-rangeR)) = circshift(GSig(rangeG,tStep+1),pkShift+glomShift);
            end
            
            % Sort the data by the rotational velocity
            if vRot(tStep) == 0
                dataRedLPB{flyID}.Stop = horzcat(dataRedLPB{flyID}.Stop, RSigNowLPB);
                dataRedRPB{flyID}.Stop = horzcat(dataRedRPB{flyID}.Stop, RSigNowRPB);
                dataGreenLPB{flyID}.Stop = horzcat(dataGreenLPB{flyID}.Stop, GSigNowLPB);
                dataGreenRPB{flyID}.Stop = horzcat(dataGreenRPB{flyID}.Stop, GSigNowRPB);
            elseif vRot(tStep) > 0
                dataRedLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRedLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    RSigNowLPB);
                dataRedRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRedRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    RSigNowRPB);
                dataGreenLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGreenLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    GSigNowLPB);
                dataGreenRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGreenRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    GSigNowRPB);
            else
                dataRedLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRedLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    RSigNowLPB);
                dataRedRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRedRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    RSigNowRPB);
                dataGreenLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGreenLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    GSigNowLPB);
                dataGreenRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGreenRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    GSigNowRPB);
            end
        end
    end
end
    
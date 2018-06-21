function [dataR, dataG] = TurnSort(allFlyData, Min, Max, Span, RRange, GRange,trial)
    
%function to sort data into bins

%input:
%allFlyData = the data to be sorted
%Min = the minimum velocity bin (deg)
%Max = the maximum velocity bin (deg)
%Span = the span of each bin (deg)
%range = the range of ROIs to be analyzed
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
            if vRot(tStep) == 0
                dataR{flyID}.Stop = horzcat(dataR{flyID}.Stop, max(RSig(RRange,tStep+1)));
                dataG{flyID}.Stop = horzcat(dataG{flyID}.Stop, max(GSig(GRange,tStep+1)));

            elseif vRot(tStep) > 0
                dataR{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataR{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    max(RSig(RRange,tStep+1)));
                dataG{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataG{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    max(GSig(GRange,tStep+1)));
            else
                dataR{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataR{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    max(RSig(RRange,tStep+1)));
                dataG{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataG{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    max(GSig(GRange,tStep+1)));
            end
        end
    end
end
    
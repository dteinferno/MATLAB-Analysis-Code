function [ROIaveMaxIntEven,ROIaveMaxIntOdd,refROIaveMaxIntOdd,ROIaveMeanEven,ROIaveMeanOdd,refROIaveMeanOdd] = ROIEvenOdd(fIDs,fIDNow,ROIs,refROIs)

info = imfinfo(fIDs(fIDNow).name);
num_images = numel(info);
recSpecs = info(1).ImageDescription;
planeLoc = strfind(recSpecs, 'numSlices');
num_planes = str2num(recSpecs(planeLoc+12));

% Read in the stack
% Get the file info and # of planes and frames
width = info(1).Width;
height = info(1).Height;
numFrames = floor(num_images/(2*num_planes));

% Find the mean volume projection for the G channel
stackMean = double(zeros(height,width,numFrames));
for incIm = 1:num_images/2
    stackMean(:,:,ceil((incIm)/num_planes)) = double(imread(fIDs(fIDNow).name, 2*incIm, 'Info', info))+ stackMean(:,:,ceil(incIm/num_planes));
end

% Find the max intensity projection
stackMaxInt = double(zeros(height,width,numFrames));
miniStack = double(zeros(height,width,num_planes));
for incIm = 1:num_images/2
    miniStack(:,:,1+mod(incIm-1,num_planes)) = double(imread(fIDs(fIDNow).name, 2*incIm, 'Info', info));
    if (1+mod(incIm-1,num_planes) == num_planes)
        stackMaxInt(:,:,ceil(incIm/num_planes)) = max(miniStack,[],3);
    end
end

stackMeanEven = zeros(size(stackMean));
stackMeanEven(2:2:end,:,:) = stackMean(2:2:end,:,:);
stackMeanOdd = zeros(size(stackMean));
stackMeanOdd(1:2:end,:,:) = stackMean(1:2:end,:,:);
stackMaxIntEven = zeros(size(stackMaxInt));
stackMaxIntEven(2:2:end,:,:) = stackMaxInt(2:2:end,:,:);
stackMaxIntOdd = zeros(size(stackMaxInt));
stackMaxIntOdd(1:2:end,:,:) = stackMean(1:2:end,:,:);

ROIaveREFMeanEven = zeros(1,size(stackMean,3));
refROIaveREFMeanEven = zeros(1,size(stackMean,3));
ROIaveREFMeanOdd = zeros(1,size(stackMean,3));
refROIaveREFMeanOdd = zeros(1,size(stackMean,3));
ROIaveREFMaxIntEven = zeros(1,size(stackMaxInt,3));
refROIaveREFMaxIntEven = zeros(1,size(stackMaxInt,3));
ROIaveREFMaxIntOdd = zeros(1,size(stackMaxInt,3));
refROIaveREFMaxIntOdd = zeros(1,size(stackMaxInt,3));

for i = 1:size(stackMean,3)
        AEven = squeeze(stackMeanEven(:,:,i));
        ROI = AEven(logical(squeeze(ROIs)));
        refROI = AEven(logical(squeeze(refROIs)));
        ROIaveREFMeanEven(i) = mean2(ROI);
        refROIaveREFMeanEven(i) = mean2(refROI);
        AOdd = squeeze(stackMeanOdd(:,:,i));
        ROI = AOdd(logical(squeeze(ROIs)));
        refROI = AOdd(logical(squeeze(refROIs)));
        ROIaveREFMeanOdd(i) = mean2(ROI);
        refROIaveREFMeanOdd(i) = mean2(refROI);
        BEven= squeeze(stackMaxIntEven(:,:,i));
        ROI = BEven(logical(squeeze(ROIs)));
        refROI = BEven(logical(squeeze(refROIs)));
        ROIaveREFMaxIntEven(i) = mean2(ROI);
        refROIaveREFMaxIntEven(i) = mean2(refROI);
        BOdd= squeeze(stackMaxIntOdd(:,:,i));
        ROI = BOdd(logical(squeeze(ROIs)));
        refROI = BOdd(logical(squeeze(refROIs)));
        ROIaveREFMaxIntOdd(i) = mean2(ROI);
        refROIaveREFMaxIntOdd(i) = mean2(refROI);
end

ROIaveMeanEven = zeros(1,length(stackMean));
ROIaveMeanOdd = zeros(1,length(stackMean));
refROIaveMeanEven = zeros(1,length(stackMean));
refROIaveMeanOdd = zeros(1,length(stackMean));

ROIaveMaxIntEven = zeros(1,length(stackMaxInt));
ROIaveMaxIntOdd = zeros(1,length(stackMaxInt));
refROIaveMaxIntEven = zeros(1,length(stackMaxInt));
refROIaveMaxIntOdd = zeros(1,length(stackMaxInt));

ROILowMeanEven = sort(squeeze(ROIaveREFMeanEven));
ROILowMeanOdd = sort(squeeze(ROIaveREFMeanOdd));
refROILowMeanEven = sort(squeeze(refROIaveREFMeanEven));
refROILowMeanOdd = sort(squeeze(refROIaveREFMeanOdd));

ROILowMaxIntEven = sort(squeeze(ROIaveREFMaxIntEven));
ROILowMaxIntOdd = sort(squeeze(ROIaveREFMaxIntOdd));
refROILowMaxIntEven = sort(squeeze(refROIaveREFMaxIntEven));
refROILowMaxIntOdd = sort(squeeze(refROIaveREFMaxIntOdd));

ROIaveMeanEven = ROIaveREFMeanEven;%./sum(ROILowMeanEven(1:floor(numFrames/10)));
ROIaveMeanOdd = ROIaveREFMeanOdd;%./sum(ROILowMeanOdd(1:floor(numFrames/10)));
refROIaveMeanEven = refROIaveREFMeanEven;%./sum(refROILowMeanEven(1:floor(numFrames/10)));
refROIaveMeanOdd = refROIaveREFMeanOdd;%./sum(refROILowMeanOdd(1:floor(numFrames/10)));

ROIaveMaxIntEven = ROIaveREFMaxIntEven;%./sum(ROILowMaxIntEven(1:floor(numFrames/10)));
ROIaveMaxIntOdd = ROIaveREFMaxIntOdd;%./sum(ROILowMaxIntOdd(1:floor(numFrames/10)));
refROIaveMaxIntEven = refROIaveREFMaxIntEven;%./sum(refROILowMaxIntEven(1:floor(numFrames/10)));
refROIaveMaxIntOdd = refROIaveREFMaxIntOdd;%./sum(refROILowMaxIntOdd(1:floor(numFrames/10)));   

end

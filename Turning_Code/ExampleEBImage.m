%% Load the stacks and the ROI data
clear;
clc;

% Specify file info
allPathname = 'D:\Imaging\2Color\20160504\';
tifName = 'Fly2_8day_6fx60D05_jRGC1ax37F06_Dark_00003_00001.tif';
moveName = 'Fly2_8day_6fx60D05_jRGC1ax37F06_Dark_03.TXT';
fileName = 'Fly2_8day_6fx60D05_jRGC1ax37F06_Dark_00003.mat';

% Load the data
tifChunks = 1;
[fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = ...
    VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,1);
load(strcat(allPathname,fileName));
load(strcat(allPathname,fileName(1:end-4),'_ROIs.mat'));

%% Gaussian Filter

% Gaussian Filter
gaussianSize = [10 10];
gaussianSigma = 4;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
GstackXYfiltMax = double(zeros(size(GstackMaxInt)));

RstackXYfiltMean = double(zeros(size(RstackMean)));
GstackXYfiltMean = double(zeros(size(GstackMean)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackMaxInt,3)
    if mod(i,100)==0
        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
    end
    RstackXYfiltMax(:,:,i) = imfilter(RstackMaxInt(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMax(:,:,i) = imfilter(GstackMaxInt(:,:,i),Gxy,'replicate','conv');
    RstackXYfiltMean(:,:,i) = imfilter(RstackMean(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMean(:,:,i) = imfilter(GstackMean(:,:,i),Gxy,'replicate','conv');
end
delete(h);

%%
ExAct = figure; 
hold on;
imagesc(fliplr(squeeze(GstackXYfiltMax(:,:,559))));
axis equal;
caxis([0 200]);
colormap('hot');
colorbar;

set(ExAct,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(ExAct,'D:\Imaging\2Color\ExAct','-dpdf');
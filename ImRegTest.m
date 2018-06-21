%% Test image registration algorithm

%% Load the data
allPathnameNow = 'C:\Users\turnerevansd\Downloads\Data\Shi\SS27\20180424\';
tifName = 'Fly1_4-5day_6fx25957_Shix27_All_00002.tif';

[stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,0);

%% Image correct the stacks
stackReg = imRegSimple(stackMaxInt, 10);

%% Image correct the stacks
stackReg2 = imRegSimple(stackReg, 10);

%% Plot the mean original stack and the mean registered stack

figure;
subplot(2,2,1);
imagesc(squeeze(mean(stackMaxInt,3)));
axis equal;
caxis([0 100]);

subplot(2,2,2);
imagesc(squeeze(mean(stackReg,3)));
axis equal;
caxis([0 100]);

subplot(2,2,3);
imagesc(squeeze(mean(stackReg2,3)));
axis equal;
caxis([0 100]);

subplot(2,2,4);
imagesc(squeeze(mean(stackReg,3))-squeeze(mean(stackReg2,3)));
axis equal;
caxis([-5 5]);

%% Make a movie of the three

% Open the movie object and make the figure
writerObj = VideoWriter(strcat(tifName(1:end-4),'.avi'));
writerObj.FrameRate= 24;
open(writerObj);
frameNow = figure('Position',[50 50 1000 1000],'Color','k');

% Generate the movie
for fm = 1:length(stackMaxInt)

    subplot(2,2,1);
    imagesc(squeeze(stackMaxInt(:,:,fm)));
    rectangle('Position',[256/2-100 256/2-100 200 200],'Curvature',1,'EdgeColor','w','LineStyle','--');
    axis equal;
    axis off;
    caxis([0 450]);

    subplot(2,2,2);
    imagesc(squeeze(stackReg(:,:,fm)));
    rectangle('Position',[256/2-100 256/2-100 200 200],'Curvature',1,'EdgeColor','w','LineStyle','--');
    axis equal;
    axis off;
    caxis([0 450]);

    subplot(2,2,3);
    imagesc(squeeze(stackReg2(:,:,fm)));
    rectangle('Position',[256/2-100 256/2-100 200 200],'Curvature',1,'EdgeColor','w','LineStyle','--');
    axis equal;
    axis off;
    caxis([0 450]);

    frame = getframe(frameNow);
    writeVideo(writerObj,frame);
end

delete(frameNow);
close(writerObj);
%% Load the stacks and the ROI data
clear;
clc;

% Specify file info
allPathname = 'D:\Imaging\2Color\20160505\';
tifName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003_00001.tif';
moveName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_03.TXT';
fileName = 'Fly2_8day_6fx37F06_jRGC1ax60D05_Dark_00003.mat';
aviName = 'Fly2_Dark_03.avi';

% Load the data
tifChunks = 1;
[fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = ...
    VRImDatLoadDiscardNoPrompt2Color(tifName,allPathname,moveName,allPathname,tifChunks,1);
load(strcat(allPathname,fileName));
load(strcat(allPathname,fileName(1:end-4),'_ROIs.mat'));

%% Gaussian Filter

% Gaussian Filter
gaussianSize = [5 5];
gaussianSigma = 2;
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

%% Make the movie - 'raw' flourescence

% Specify time periods to consider for making the movie
movieTimes = zeros(4,2);
movieTimes(1,1) = 6.25;
movieTimes(1,2) = 18.25;
movieTimes(2,1) = 21;
movieTimes(2,2) = 31;
movieTimes(3,1) = 41;
movieTimes(3,2) = 51;
movieTimes(4,1) = 56;
movieTimes(4,2) = 66;

% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GstackXYfiltMax,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);
vRotNorm = vRot./max(abs(vRot));

% Savitzky-Golay Filter the RFs
sgolayOrder = 3;
sgolayWindow = 11;
RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);

% Find the PVA
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';
for ts = 1:length(GROIaveMaxFilt)
    meanAngRawG(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanAngRawR(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    meanIntRawG(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanIntRawR(ts) = circ_r(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
end

for tFrame = 1:size(movieTimes,1)
    
    % Find the movie start and stop times
    tStart = movieTimes(tFrame,1);
    tStop = movieTimes(tFrame,2);
    minFG = max(minFGAll,round(tStart*framerate));
    maxFG = min(maxFGAll,round(tStop*framerate));
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(moveName(1:end-4),'_',num2str(tStart),'to',num2str(tStop),'.avi'));
    writerObj.FrameRate= framerate;
    open(writerObj);

    % Create the movie figure
    movieFig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','k');

    % Set the colorbar limits
    RmaxVal = max(max(max(RstackXYfiltMax(:,:,minFG:maxFG))));
    RminVal = RmaxVal/10;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GstackXYfiltMax(:,:,minFG:maxFG))));
    GminVal = GmaxVal/20;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    for clrNow = 1:3
    
        % Get the current frame
        curRFrame = (RstackXYfiltMax(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
        curGFrame = 3*(GstackXYfiltMax(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

        % Create the colored image based off of the frame
        RIm = zeros([size(curRFrame) 3]);
        RIm(:,:,1) = curRFrame;
        RIm(:,:,3) = curRFrame;
        GIm = zeros([size(curGFrame) 3]);
        GIm(:,:,2) = curGFrame;
        overlayIm = zeros([size(curRFrame) 3]);
        overlayIm(:,:,1) = curRFrame;
        overlayIm(:,:,2) = curGFrame;
        overlayIm(:,:,3) = curRFrame;
        
        if clrNow == 1
            imageNow = RIm;
        elseif clrNow == 2
            imageNow = GIm;
        else
            imageNow = overlayIm;
        end
        imageNow = flipud(fliplr(imageNow));

        % Create the base figure;
        subplot(6,4,[1:16]);
        hold on;
        imshow(imageNow);
        if clrNow == 1
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
            line([25 50], [252 252],'Color','m','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
        elseif clrNow == 2
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([25 50], [252 252],'Color','g','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
        else
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
        end
        axis off;
        if clrNow == 1
            title('E-PG fluorescence','Color','w','FontSize',14);
        elseif clrNow == 2
            title('P-EN fluorescence','Color','w','FontSize',14);
        else
            title('Combined fluorescence','Color','w','FontSize',14);
        end

        % Look at the fly on the ball
        subplot(6,4,[17 18 21 22]);
        walkVid = VideoReader(aviName);
        vidWidth = walkVid.Width;
        vidHeight = walkVid.Height;
        mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
            'colormap',[]);
        mov(1).cdata = readFrame(walkVid);
        for advmovie = 2:floor(minFG*num_planes/2)-1
            mov(advmovie).cdata = readFrame(walkVid);
        end
        imshow(mov(advmovie).cdata(50:300,100:385,1));
        
        % Track heading or rotational velocity
        subplot(6,4,[19 20 23 24]);
        cla;
        if clrNow == 1 || clrNow == 3
            subplot(6,4,[19 20 23 24]);
            cla;
            hold on;
            rectangle('Position',[0 0 256 256],'EdgeColor','w','Curvature',1);
            line([128 128+128*cos(OffsetRotMatch(minFG-minFGAll+1,2))], [128 128-128*sin(OffsetRotMatch(minFG-minFGAll+1,2))],'Color','w','LineWidth',2);
            axis off;
            axis equal;
            title('heading','Color','w','FontSize',14);
        else
            subplot(6,4,[19 23]);
            cla;
            bar(abs(vRotNorm(minFG-minFGAll+1)),'FaceColor','w');
            ylim([0 1]);
            axis off;
            title('|v_{rot}|','Color','w','FontSize',14);
        end
        
        frame = getframe(movieFig);
        writeVideo(writerObj,frame);

        % Hold steady
        tHold = 3;
        for frameHold = 1:round(tHold*framerate)
            writeVideo(writerObj,frame);
        end

        for frameNum = minFG+1:maxFG-1

            curRFrame = (RstackXYfiltMax(:,:,frameNum)-RminVal)./(RmaxVal-RminVal);
            curGFrame = 3*(GstackXYfiltMax(:,:,frameNum)-GminVal)./(GmaxVal-GminVal);

            RIm = zeros([size(curRFrame) 3]);
            RIm(:,:,1) = curRFrame;
            RIm(:,:,3) = curRFrame;
            GIm = zeros([size(curGFrame) 3]);
            GIm(:,:,2) = curGFrame;
            overlayIm = zeros([size(curRFrame) 3]);
            overlayIm(:,:,1) = curRFrame;
            overlayIm(:,:,2) = curGFrame;
            overlayIm(:,:,3) = curRFrame;
            
            if clrNow == 1
                imageNow = RIm;
            elseif clrNow == 2
                imageNow = GIm;
            else
                imageNow = overlayIm;
            end
            imageNow = flipud(fliplr(imageNow));

           % Create the base figure;
            subplot(6,4,[1:16]);
            cla;
            hold on;
            imshow(imageNow);
            if clrNow == 1
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
            elseif clrNow == 2
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
            else
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
            end
            
            axis off;
            text(200, 240, strcat('t=',num2str(OffsetRotMatch(frameNum-minFGAll,1)-OffsetRotMatch(minFG+1-minFGAll,1),3),'s'),'Color','w','FontSize',14);

            % Look at the fly on the ball            
            subplot(6,4,[17 18 21 22]);
            mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
                'colormap',[]);
            for advmovie = floor((frameNum-1)*num_planes/2):floor(frameNum*num_planes/2)-1
                mov(advmovie).cdata = readFrame(walkVid);
            end
            imshow(mov(advmovie).cdata(50:300,100:385,1));

            % Track heading or rotational velocity
            if clrNow == 1 || clrNow == 3
                subplot(6,4,[19 20 23 24]);
                cla;
                hold on;
                rectangle('Position',[0 0 256 256],'EdgeColor','w','Curvature',1);
                line([128 128+128*cos(OffsetRotMatch(frameNum-minFGAll+1,2))], [128 128-128*sin(OffsetRotMatch(frameNum-minFGAll+1,2))],'Color','w','LineWidth',2);
                axis off;
                axis equal;
                title('heading','Color','w','FontSize',14);
            else
                subplot(6,4,[19 23]);
                cla;
                bar(abs(vRotNorm(frameNum-minFGAll+1)),'FaceColor','w');
                ylim([0 1]);
                axis off;
                title('|v_{rot}|','Color','w','FontSize',14);
            end
                        
            frame = getframe(movieFig);
            writeVideo(writerObj,frame);
        end
    end

    delete(movieFig);
    close(writerObj);
end

%% Make the movie - background subtracted flourescence
RStackSub = zeros(size(RstackMaxInt));
GStackSub = zeros(size(GstackMaxInt));
RMean = mean(RstackMaxInt,3);
GMean = mean(GstackMaxInt,3);
for frame = 1:size(RstackMaxInt,3)
    RStackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
    GStackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
end

% Gaussian Filter
gaussianSize = [7 7];
gaussianSigma = 3;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RStackSubFilt = double(zeros(size(RstackSub)));
GStackSubFilt = double(zeros(size(GstackSub)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackSub,3)
    if mod(i,100)==0
        waitbar(i/length(RstackSub),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackSub))]);
    end
    RStackSubFilt(:,:,i) = imfilter(RstackSub(:,:,i),Gxy,'replicate','conv');
    GStackSubFilt(:,:,i) = imfilter(GstackSub(:,:,i),Gxy,'replicate','conv');
end
delete(h);

% Specify time periods to consider for making the movie
movieTimes = zeros(4,2);
movieTimes(1,1) = 6.5;
movieTimes(1,2) = 18;
movieTimes(2,1) = 21;
movieTimes(2,2) = 31;
movieTimes(3,1) = 41;
movieTimes(3,2) = 51;
movieTimes(4,1) = 56;
movieTimes(4,2) = 66;

% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GstackXYfiltMax,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);
vRotNorm = vRot./max(abs(vRot));

% Savitzky-Golay Filter the RFs
sgolayOrder = 3;
sgolayWindow = 11;
RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);

% Find the PVA
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';
for ts = 1:length(GROIaveMaxFilt)
    meanAngRawG(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanAngRawR(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    meanIntRawG(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanIntRawR(ts) = circ_r(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
end

for tFrame = 1:size(movieTimes,1)
    
    % Find the movie start and stop times
    tStart = movieTimes(tFrame,1);
    tStop = movieTimes(tFrame,2);
    minFG = max(minFGAll,round(tStart*framerate));
    maxFG = min(maxFGAll,round(tStop*framerate));
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(moveName(1:end-4),'_',num2str(tStart),'to',num2str(tStop),'_subtracted.avi'));
    writerObj.FrameRate= framerate;
    open(writerObj);

    % Create the movie figure
    movieFig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','k');

    % Set the colorbar limits
    RmaxVal = max(max(max(RStackSubFilt(:,:,minFG:maxFG))));
    RminVal = RmaxVal/10;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GStackSubFilt(:,:,minFG:maxFG))));
    GminVal = GmaxVal/20;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    for clrNow = 1:3
    
        % Get the current frame
        curRFrame = 2*(RStackSubFilt(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
        curGFrame = 3*(GStackSubFilt(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

        % Create the colored image based off of the frame
        RIm = zeros([size(curRFrame) 3]);
        RIm(:,:,1) = curRFrame;
        RIm(:,:,3) = curRFrame;
        GIm = zeros([size(curGFrame) 3]);
        GIm(:,:,2) = curGFrame;
        overlayIm = zeros([size(curRFrame) 3]);
        overlayIm(:,:,1) = curRFrame;
        overlayIm(:,:,2) = curGFrame;
        overlayIm(:,:,3) = curRFrame;
        
        if clrNow == 1
            imageNow = RIm;
        elseif clrNow == 2
            imageNow = GIm;
        else
            imageNow = overlayIm;
        end
        imageNow = flipud(fliplr(imageNow));

        % Create the base figure;
        subplot(6,4,[1:16]);
        hold on;
        imshow(imageNow);
        if clrNow == 1
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
            line([25 50], [252 252],'Color','m','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
            rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[110 115 256-220 256-230],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
        elseif clrNow == 2
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([25 50], [252 252],'Color','g','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
            rectangle('Position',[40 70 256-80 256-140],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
        else
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
            rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
        end
        axis off;
        if clrNow == 1
            title('E-PG fluorescence above background','Color','w','FontSize',14);
        elseif clrNow == 2
            title('P-EN fluorescence above background','Color','w','FontSize',14);
        else
            title('Combined fluorescence above background','Color','w','FontSize',14);
        end

        % Look at the fly on the ball
        subplot(6,4,[17 18 21 22]);
        walkVid = VideoReader(aviName);
        vidWidth = walkVid.Width;
        vidHeight = walkVid.Height;
        mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
            'colormap',[]);
        mov(1).cdata = readFrame(walkVid);
        for advmovie = 2:floor(minFG*num_planes/2)-1
            mov(advmovie).cdata = readFrame(walkVid);
        end
        imshow(mov(advmovie).cdata(50:300,100:385,1));
        
        % Track heading or rotational velocity
        subplot(6,4,[19 20 23 24]);
        cla;
        if clrNow == 1 || clrNow == 3
            subplot(6,4,[19 20 23 24]);
            cla;
            hold on;
            rectangle('Position',[0 0 256 256],'EdgeColor','w','Curvature',1);
            line([128 128+128*cos(OffsetRotMatch(minFG-minFGAll+1,2))], [128 128-128*sin(OffsetRotMatch(minFG-minFGAll+1,2))],'Color','w','LineWidth',2);
            axis off;
            axis equal;
            title('heading','Color','w','FontSize',14);
        else
            subplot(6,4,[19 23]);
            cla;
            bar(abs(vRotNorm(minFG-minFGAll+1)),'FaceColor','w');
            ylim([0 1]);
            axis off;
            title('|v_{rot}|','Color','w','FontSize',14);
        end
        
        frame = getframe(movieFig);
        writeVideo(writerObj,frame);

        % Hold steady
        tHold = 3;
        for frameHold = 1:round(tHold*framerate)
            writeVideo(writerObj,frame);
        end

        for frameNum = minFG+1:maxFG-1

            curRFrame = 2*(RStackSubFilt(:,:,frameNum)-RminVal)./(RmaxVal-RminVal);
            curGFrame = 3*(GStackSubFilt(:,:,frameNum)-GminVal)./(GmaxVal-GminVal);

            RIm = zeros([size(curRFrame) 3]);
            RIm(:,:,1) = curRFrame;
            RIm(:,:,3) = curRFrame;
            GIm = zeros([size(curGFrame) 3]);
            GIm(:,:,2) = curGFrame;
            overlayIm = zeros([size(curRFrame) 3]);
            overlayIm(:,:,1) = curRFrame;
            overlayIm(:,:,2) = curGFrame;
            overlayIm(:,:,3) = curRFrame;
            
            if clrNow == 1
                imageNow = RIm;
            elseif clrNow == 2
                imageNow = GIm;
            else
                imageNow = overlayIm;
            end
            imageNow = flipud(fliplr(imageNow));

           % Create the base figure;
            subplot(6,4,[1:16]);
            cla;
            hold on;
            imshow(imageNow);
            if clrNow == 1
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
                rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[110 115 256-220 256-230],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            elseif clrNow == 2
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
                rectangle('Position',[40 70 256-80 256-140],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            else
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
                rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            end
            
            axis off;
            text(200, 240, strcat('t=',num2str(OffsetRotMatch(frameNum-minFGAll,1)-OffsetRotMatch(minFG+1-minFGAll,1),3),'s'),'Color','w','FontSize',14);

            % Look at the fly on the ball            
            subplot(6,4,[17 18 21 22]);
            mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
                'colormap',[]);
            for advmovie = floor((frameNum-1)*num_planes/2):floor(frameNum*num_planes/2)-1
                mov(advmovie).cdata = readFrame(walkVid);
            end
            imshow(mov(advmovie).cdata(50:300,100:385,1));

            % Track heading or rotational velocity
            if clrNow == 1 || clrNow == 3
                subplot(6,4,[19 20 23 24]);
                cla;
                hold on;
                rectangle('Position',[0 0 256 256],'EdgeColor','w','Curvature',1);
                line([128 128+128*cos(OffsetRotMatch(frameNum-minFGAll+1,2))], [128 128-128*sin(OffsetRotMatch(frameNum-minFGAll+1,2))],'Color','w','LineWidth',2);
                axis off;
                axis equal;
                title('heading','Color','w','FontSize',14);
            else
                subplot(6,4,[19 23]);
                cla;
                bar(abs(vRotNorm(frameNum-minFGAll+1)),'FaceColor','w');
                ylim([0 1]);
                axis off;
                title('|v_{rot}|','Color','w','FontSize',14);
            end
                        
            frame = getframe(movieFig);
            writeVideo(writerObj,frame);
        end
    end

    delete(movieFig);
    close(writerObj);
end

%% Make the movie - background subtracted flourescence w/Heading overlay
RStackSub = zeros(size(RstackMaxInt));
GStackSub = zeros(size(GstackMaxInt));
RMean = mean(RstackMaxInt,3);
GMean = mean(GstackMaxInt,3);
for frame = 1:size(RstackMaxInt,3)
    RStackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
    GStackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
end

% Gaussian Filter
gaussianSize = [7 7];
gaussianSigma = 3;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RStackSubFilt = double(zeros(size(RStackSub)));
GStackSubFilt = double(zeros(size(GStackSub)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RStackSub,3)
    if mod(i,100)==0
        waitbar(i/length(RStackSub),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RStackSub))]);
    end
    RStackSubFilt(:,:,i) = imfilter(RStackSub(:,:,i),Gxy,'replicate','conv');
    GStackSubFilt(:,:,i) = imfilter(GStackSub(:,:,i),Gxy,'replicate','conv');
end
delete(h);

% Specify time periods to consider for making the movie
movieTimes = zeros(4,2);
movieTimes(1,1) = 6.5;
movieTimes(1,2) = 11;
movieTimes(2,1) = 21;
movieTimes(2,2) = 31;
movieTimes(3,1) = 45;
movieTimes(3,2) = 51;
movieTimes(4,1) = 60;
movieTimes(4,2) = 66;

% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GStackSubFilt,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);

sgolayOrder = 3;
sgolayWindow = 11;
sgolayfilt(vRot,sgolayOrder,sgolayWindow);
vRotNorm = vRot./max(abs(vRot));

% Savitzky-Golay Filter the RFs
RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);

% Find the PVA
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';
for ts = 1:length(GROIaveMaxFilt)
    meanAngRawG(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanAngRawR(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    meanIntRawG(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanIntRawR(ts) = circ_r(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
end

for tFrame = 1:size(movieTimes,1)
    
    % Find the movie start and stop times
    tStart = movieTimes(tFrame,1);
    tStop = movieTimes(tFrame,2);
    minFG = max(minFGAll,round(tStart*framerate));
    maxFG = min(maxFGAll,round(tStop*framerate));
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(moveName(1:end-4),'_',num2str(tStart),'to',num2str(tStop),'_subtracted.avi'));
    writerObj.FrameRate= framerate;
    open(writerObj);

    % Create the movie figure
    movieFig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','k');

    % Set the colorbar limits
    RmaxVal = max(max(max(RStackSubFilt(:,:,minFG:maxFG))));
    RminVal = RmaxVal/10;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GStackSubFilt(:,:,minFG:maxFG))));
    GminVal = GmaxVal/20;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    for clrNow = 1:3
    
        % Get the current frame
        curRFrame = 2*(RStackSubFilt(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
        curGFrame = 3*(GStackSubFilt(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

        % Create the colored image based off of the frame
        RIm = zeros([size(curRFrame) 3]);
        RIm(:,:,1) = curRFrame;
        RIm(:,:,3) = curRFrame;
        GIm = zeros([size(curGFrame) 3]);
        GIm(:,:,2) = curGFrame;
        overlayIm = zeros([size(curRFrame) 3]);
        overlayIm(:,:,1) = curRFrame;
        overlayIm(:,:,2) = curGFrame;
        overlayIm(:,:,3) = curRFrame;
        
        if clrNow == 1
            imageNow = RIm;
        elseif clrNow == 2
            imageNow = GIm;
        else
            imageNow = overlayIm;
        end
        imageNow = flipud(fliplr(imageNow));

        % Create the base figure;
        subplot(6,4,[1:16]);
        cla;
        hold on;
        imshow(imageNow);
        if clrNow == 1
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
            line([25 50], [252 252],'Color','m','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
            
            rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[110 115 256-220 256-230],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            
            OffsetAng = meanAngRawR(minFG);
            rectangle('Position',[...
                    128+0.5*(256-20)*cos(OffsetAng)-5 ...
                    128+0.5*(256-110)*sin(OffsetAng)-5 ...
                    10 10],...
                    'Curvature',[1 1],'FaceColor','w');
            
            rectangle('Position',[35 230 10 10],...
                    'Curvature',[1 1],'FaceColor','w');
            text(60, 235, '= heading','Color','w','FontSize', 18);
                
        elseif clrNow == 2
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([25 50], [252 252],'Color','g','LineWidth',2);
            text(60, 250, '= population vector average','Color','w','FontSize', 18);
            
            rectangle('Position',[40 70 256-80 256-140],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            
            hvR = abs(vRotNorm(minFG-minFGAll+1));
            rectangle('Position',[245 256-70-60*hvR 15 60*hvR],'FaceColor','w');
            text(240, 256-130, '|v_{Rot}|','FontSize',18,'Color','w');
        else
            line([128 128+500*meanIntRawG(minFG)*cos(meanAngRawG(minFG))], [128 128+500*meanIntRawG(minFG)*sin(meanAngRawG(minFG))],'Color','g','LineWidth',2);
            line([128 128+500*meanIntRawR(minFG)*cos(meanAngRawR(minFG))], [128 128+500*meanIntRawR(minFG)*sin(meanAngRawR(minFG))],'Color','m','LineWidth',2);
            
            rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
            
            OffsetAng = meanAngRawR(minFG);
            rectangle('Position',[...
                    128+0.5*(256-20)*cos(OffsetAng)-5 ...
                    128+0.5*(256-110)*sin(OffsetAng)-5 ...
                    10 10],...
                    'Curvature',[1 1],'FaceColor','w');
        end
        axis off;
        if clrNow == 1
            title('E-PG fluorescence above background','Color','w','FontSize',14);
        elseif clrNow == 2
            title('P-EN fluorescence above background','Color','w','FontSize',14);
        else
            title('Combined fluorescence above background','Color','w','FontSize',14);
        end

        % Look at the fly on the ball
        subplot(6,4,[18 19 22 23]);
        walkVid = VideoReader(aviName);
        vidWidth = walkVid.Width;
        vidHeight = walkVid.Height;
        mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
            'colormap',[]);
        mov(1).cdata = readFrame(walkVid);
        for advmovie = 2:floor(minFG*num_planes/2)-1
            mov(advmovie).cdata = readFrame(walkVid);
        end
        imshow(mov(advmovie).cdata(50:300,100:385,1));
        
        frame = getframe(movieFig);
        writeVideo(writerObj,frame);

        % Hold steady
        tHold = 3;
        for frameHold = 1:round(tHold*framerate)
            writeVideo(writerObj,frame);
        end

        for frameNum = minFG+1:maxFG-1

            curRFrame = 2*(RStackSubFilt(:,:,frameNum)-RminVal)./(RmaxVal-RminVal);
            curGFrame = 3*(GStackSubFilt(:,:,frameNum)-GminVal)./(GmaxVal-GminVal);

            RIm = zeros([size(curRFrame) 3]);
            RIm(:,:,1) = curRFrame;
            RIm(:,:,3) = curRFrame;
            GIm = zeros([size(curGFrame) 3]);
            GIm(:,:,2) = curGFrame;
            overlayIm = zeros([size(curRFrame) 3]);
            overlayIm(:,:,1) = curRFrame;
            overlayIm(:,:,2) = curGFrame;
            overlayIm(:,:,3) = curRFrame;
            
            if clrNow == 1
                imageNow = RIm;
            elseif clrNow == 2
                imageNow = GIm;
            else
                imageNow = overlayIm;
            end
            imageNow = flipud(fliplr(imageNow));

           % Create the base figure;
            subplot(6,4,[1:16]);
            cla;
            hold on;
            imshow(imageNow);
            if clrNow == 1
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
                
                rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[110 115 256-220 256-230],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                
                OffsetAng = OffsetRotMatch(frameNum-minFGAll+1,2) - OffsetRotMatch(minFG-minFGAll+1,2) + meanAngRawR(minFG);
                rectangle('Position',[...
                    128+0.5*(256-20)*cos(OffsetAng) - 5 ...
                    128+0.5*(256-110)*sin(OffsetAng) - 5 ...
                    10 10],...
                    'Curvature',[1 1],'FaceColor','w');
            elseif clrNow == 2
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
                rectangle('Position',[40 70 256-80 256-140],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                
                hvR = abs(vRotNorm(frameNum-minFGAll+1));
                rectangle('Position',[245 256-70-60*hvR 15 60*hvR],'FaceColor','w');
                text(240, 256-130, '|v_{Rot}|','FontSize',18,'Color','w');
            else
                line([128 128+500*meanIntRawG(frameNum)*cos(meanAngRawG(frameNum))], [128 128+500*meanIntRawG(frameNum)*sin(meanAngRawG(frameNum))],'Color','g','LineWidth',2);
                line([128 128+500*meanIntRawR(frameNum)*cos(meanAngRawR(frameNum))], [128 128+500*meanIntRawR(frameNum)*sin(meanAngRawR(frameNum))],'Color','m','LineWidth',2);
                
                rectangle('Position',[10 55 256-20 256-110],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                rectangle('Position',[120 120 256-240 256-240],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5,'LineStyle','--');
                
                OffsetAng = OffsetRotMatch(frameNum-minFGAll+1,2) - OffsetRotMatch(minFG-minFGAll+1,2) + meanAngRawR(minFG);
                rectangle('Position',[...
                    128+0.5*(256-20)*cos(OffsetAng) - 5 ...
                    128+0.5*(256-110)*sin(OffsetAng) - 5 ...
                    10 10],...
                    'Curvature',[1 1],'FaceColor','w');
            end
            
            axis off;
            text(200, 240, strcat('t=',num2str(OffsetRotMatch(frameNum-minFGAll,1)-OffsetRotMatch(minFG+1-minFGAll,1),3),'s'),'Color','w','FontSize',14);

            % Look at the fly on the ball            
            subplot(6,4,[18 19 22 23]);
            mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
                'colormap',[]);
            for advmovie = floor((frameNum-1)*num_planes/2):floor(frameNum*num_planes/2)-1
                mov(advmovie).cdata = readFrame(walkVid);
            end
            imshow(mov(advmovie).cdata(50:300,100:385,1));
                        
            frame = getframe(movieFig);
            writeVideo(writerObj,frame);
        end
    end

    delete(movieFig);
    close(writerObj);
end

%% Subtract the background

RStackSub = zeros(size(RstackMaxInt));
GStackSub = zeros(size(GstackMaxInt));
RMean = mean(RstackMaxInt,3);
GMean = mean(GstackMaxInt,3);
for frame = 1:size(RstackMaxInt,3)
    RStackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
    GStackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
end

% Gaussian Filter
gaussianSize = [25 25];
gaussianSigma = 10;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

RStackSubFilt = double(zeros(size(RStackSub)));
GStackSubFilt = double(zeros(size(GStackSub)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RStackSub,3)
    if mod(i,100)==0
        waitbar(i/length(RStackSub),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RStackSub))]);
    end
    RStackSubFilt(:,:,i) = imfilter(RStackSub(:,:,i),Gxy,'replicate','conv');
    GStackSubFilt(:,:,i) = imfilter(GStackSub(:,:,i),Gxy,'replicate','conv');
end
delete(h);

%% Plot three frames of interest
tPts(1) = 61.35;
tPts(2) = 61.5;
tPts(3) = 61.8;

% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GStackSubFilt,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);
vRotNorm = vRot./max(abs(vRot));

PrettyPic = figure;

for tNow = 1:length(tPts)

    minFG = max(minFGAll,round(tPts(tNow)*framerate));
    maxFG = min(maxFGAll,round(tPts(3)*framerate));
    
    % Set the colorbar limits
    RmaxVal = max(max(max(RStackSubFilt(:,:,minFG:maxFG))));
    RminVal = RmaxVal/5;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GStackSubFilt(:,:,minFG:maxFG))));
    GminVal = GmaxVal/12;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    % Get the current frame
    curRFrame = 5*(RStackSubFilt(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
    curGFrame = 8*(GStackSubFilt(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

    % Create the colored image based off of the frame
    RIm = zeros([size(curRFrame) 3]);
    RIm(:,:,1) = curRFrame;
    RIm(:,:,3) = curRFrame;
    GIm = zeros([size(curGFrame) 3]);
    GIm(:,:,2) = curGFrame;
    overlayIm = zeros([size(curRFrame) 3]);
    overlayIm(:,:,1) = curRFrame;
    overlayIm(:,:,2) = curGFrame;
    overlayIm(:,:,3) = curRFrame;

    imageNow = zeros(2*size(RIm,1),2*size(RIm,1),size(RIm,3));
    if tNow == 1
        imageNow(1:2:end,1:2:end,:) = RIm;
    elseif tNow == 2
        imageNow(1:2:end,1:2:end,:) = overlayIm;
    else
        imageNow(1:2:end,1:2:end,:) = RIm;
    end
    
    sz = 4;
    H = fspecial('average',[sz sz]);
    imageNow(:,:,1) = conv2(squeeze(imageNow(:,:,1)),H,'same');
    imageNow(:,:,2) = conv2(squeeze(imageNow(:,:,2)),H,'same');
    imageNow(:,:,3) = conv2(squeeze(imageNow(:,:,3)),H,'same');
    
    imageNow = flipud(fliplr(imageNow));

    % Create the base figure;
    subplot(1,3,tNow);
    hold on;
    imshow(imageNow);
    rectangle('Position',[80 100 512-160 512-200],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5);
    rectangle('Position',[230 230 512-460 512-460],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5);
    xlim([50 462]);
    ylim([50 462]);
    
end

set(PrettyPic,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PrettyPic,'D:\Imaging\2Color\PrettyPic','-dpdf');

%% Make the movie - background subtracted flourescence - for Yi

% Specify time periods to consider for making the movie
movieTimes = zeros(4,2);
movieTimes(1,1) = 6;
movieTimes(1,2) = 66;

% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GstackXYfiltMax,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);

% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);
vRotNorm = vRot./max(abs(vRot));

% Savitzky-Golay Filter the RFs
sgolayOrder = 3;
sgolayWindow = 11;
RROIaveMaxFilt = sgolayfilt(RROIaveMax,sgolayOrder,sgolayWindow,[],2);
GROIaveMaxFilt = sgolayfilt(GROIaveMax,sgolayOrder,sgolayWindow,[],2);

% Find the PVA
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';
for ts = 1:length(GROIaveMaxFilt)
    meanAngRawG(ts) = circ_mean(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanAngRawR(ts) = circ_mean(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
    meanIntRawG(ts) = circ_r(angsraw, squeeze(GROIaveMaxFilt(:,ts)));
    meanIntRawR(ts) = circ_r(angsraw, squeeze(RROIaveMaxFilt(:,ts)));
end

for tFrame = 1:size(movieTimes,1)
    
    % Find the movie start and stop times
    tStart = movieTimes(tFrame,1);
    tStop = movieTimes(tFrame,2);
    minFG = max(minFGAll,round(tStart*framerate));
    maxFG = min(maxFGAll,round(tStop*framerate));
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(moveName(1:end-4),'_',num2str(tStart),'to',num2str(tStop),'_subtracted_forYi.avi'));
    writerObj.FrameRate= framerate;
    open(writerObj);

    % Create the movie figure
    movieFig = figure('units','normalized','outerposition',[0 0 0.5 1],'Color','k');

    % Set the colorbar limits
    RmaxVal = max(max(max(RStackSubFilt(:,:,minFG:maxFG))));
    RminVal = RmaxVal/10;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GStackSubFilt(:,:,minFG:maxFG))));
    GminVal = GmaxVal/20;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    for clrNow = 3
    
        % Get the current frame
        curRFrame = 2*(RStackSubFilt(:,:,minFG)-RminVal)./(RmaxVal-RminVal);
        curGFrame = 3*(GStackSubFilt(:,:,minFG)-GminVal)./(GmaxVal-GminVal);

        % Create the colored image based off of the frame
        RIm = zeros([size(curRFrame) 3]);
        RIm(:,:,1) = curRFrame;
        RIm(:,:,3) = curRFrame;
        GIm = zeros([size(curGFrame) 3]);
        GIm(:,:,2) = curGFrame;
        overlayIm = zeros([size(curRFrame) 3]);
        overlayIm(:,:,1) = curRFrame;
        overlayIm(:,:,2) = curGFrame;
        overlayIm(:,:,3) = curRFrame;
        
        if clrNow == 1
            imageNow = RIm;
        elseif clrNow == 2
            imageNow = GIm;
        else
            imageNow = overlayIm;
        end
        imageNow = flipud(fliplr(imageNow));

        % Create the base figure;
        subplot(4,4,[1:16]);
        hold on;
        imshow(imageNow);
        axis off;
        text(200, 240, strcat('t=',num2str(OffsetRotMatch(minFG-minFGAll,1)-OffsetRotMatch(minFG+1-minFGAll,1),3),'s'),'Color','w','FontSize',14);
        text(35, 230, '20\mum','Color','w','FontSize',10);
        line([20 966],[245 245],'Color','w','LineWidth',2);
        
        frame = getframe(movieFig);
        writeVideo(writerObj,frame);

        for frameNum = minFG+1:maxFG-1

            curRFrame = 2*(RStackSubFilt(:,:,frameNum)-RminVal)./(RmaxVal-RminVal);
            curGFrame = 3*(GStackSubFilt(:,:,frameNum)-GminVal)./(GmaxVal-GminVal);

            RIm = zeros([size(curRFrame) 3]);
            RIm(:,:,1) = curRFrame;
            RIm(:,:,3) = curRFrame;
            GIm = zeros([size(curGFrame) 3]);
            GIm(:,:,2) = curGFrame;
            overlayIm = zeros([size(curRFrame) 3]);
            overlayIm(:,:,1) = curRFrame;
            overlayIm(:,:,2) = curGFrame;
            overlayIm(:,:,3) = curRFrame;
            
            if clrNow == 1
                imageNow = RIm;
            elseif clrNow == 2
                imageNow = GIm;
            else
                imageNow = overlayIm;
            end
            imageNow = flipud(fliplr(imageNow));

           % Create the base figure;
            subplot(4,4,[1:16]);
            cla;
            hold on;
            imshow(imageNow);
            axis off;
            text(200, 240, strcat('t=',num2str(OffsetRotMatch(frameNum-minFGAll,1)-OffsetRotMatch(minFG+1-minFGAll,1),3),'s'),'Color','w','FontSize',14);
            text(35, 230, '20\mum','Color','w','FontSize',10);
            line([20 96],[245 245],'Color','w','LineWidth',2);         
            frame = getframe(movieFig);
            writeVideo(writerObj,frame);
        end
    end

    delete(movieFig);
    close(writerObj);
end

%% Make small movies for publicity
tStart = 61.2;
tStop = 62.3;
    
% Find the number of planes and the framerate
num_planes = length(positionDat.tFrameGrab)./size(GStackSubFilt,3);
framerate = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
minFG = max(minFGAll,round(tStart*framerate));
maxFG = min(maxFGAll,round(tStop*framerate));
minFGAll = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);  
maxFGAll = round(length(positionDat.tFrameGrab)/num_planes);


% Match the position data to the framegrab times
OffsetRotMatch = zeros(maxFGAll-minFGAll+1,2);
OffsetForMatch = zeros(maxFGAll-minFGAll+1,1);
OffsetLatMatch = zeros(maxFGAll-minFGAll+1,1);
for interp = minFGAll:maxFGAll
    tMatch = find(positionDat.t >= (positionDat.t(1) +...
        (positionDat.tFrameGrab((interp-1)*num_planes+1)- positionDat.tVR(1))/10000));
    OffsetRotMatch(interp-minFGAll+1,1) = positionDat.t(tMatch(1));
    OffsetRotMatch(interp-minFGAll+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
    OffsetForMatch(interp-minFGAll+1) = positionDat.OffsetFor(tMatch(1));
    OffsetLatMatch(interp-minFGAll+1) = positionDat.OffsetLat(tMatch(1));
end
OffsetRotMatchUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
vRot = diff(OffsetRotMatchUnwrap);
vRotNorm = vRot./max(abs(vRot));

movieFig = figure;
% Open the movie object and make the figure
writerObj = VideoWriter('Lturn_R.avi');
writerObj.FrameRate= framerate;
open(writerObj);

for tNow = 1:maxFG-minFG-1

    % Set the colorbar limits
    RmaxVal = max(max(max(RStackSubFilt(:,:,minFG:maxFG))));
    RminVal = RmaxVal/5;%4*min(min(min(RstackXYfiltMax(:,:,minFG:maxFG))));
    GmaxVal = max(max(max(GStackSubFilt(:,:,minFG:maxFG))));
    GminVal = GmaxVal/12;%2*min(min(min(GstackXYfiltMax(:,:,minFG:maxFG))));
        
    % Get the current frame
    curRFrame = 5*(RStackSubFilt(:,:,tNow+minFG-1)-RminVal)./(RmaxVal-RminVal);
    curGFrame = 8*(GStackSubFilt(:,:,tNow+minFG-1)-GminVal)./(GmaxVal-GminVal);

    % Create the colored image based off of the frame
    RIm = zeros([size(curRFrame) 3]);
    RIm(:,:,1) = curRFrame;
    RIm(:,:,3) = curRFrame;
    GIm = zeros([size(curGFrame) 3]);
    GIm(:,:,2) = curGFrame;
    overlayIm = zeros([size(curRFrame) 3]);
    overlayIm(:,:,1) = curRFrame;
    overlayIm(:,:,2) = curGFrame;
    overlayIm(:,:,3) = curRFrame;

    for clr = 1
        imageNow = zeros(2*size(RIm,1),2*size(RIm,1),size(RIm,3));
        if clr == 1
            imageNow(1:2:end,1:2:end,:) = RIm;
        else
            imageNow(1:2:end,1:2:end,:) = overlayIm;
        end

        sz = 4;
        H = fspecial('average',[sz sz]);
        imageNow(:,:,1) = conv2(squeeze(imageNow(:,:,1)),H,'same');
        imageNow(:,:,2) = conv2(squeeze(imageNow(:,:,2)),H,'same');
        imageNow(:,:,3) = conv2(squeeze(imageNow(:,:,3)),H,'same');

        imageNow = flipud(fliplr(imageNow));

        % Create the base figure;
%         subplot(2,maxFG-minFG-1,tNow+(clr-1)*(maxFG-minFG-1));
        hold on;
        imshow(imageNow);
        rectangle('Position',[80 100 512-160 512-200],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5);
        rectangle('Position',[230 230 512-460 512-460],'Curvature',[1 1],'EdgeColor','w','LineWidth',0.5);
        xlim([50 462]);
        ylim([50 462]);
    end
    
    frame = getframe(movieFig);
    writeVideo(writerObj,frame);
    
end

 delete(movieFig);
    close(writerObj);

% set(PrettyPic,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(PrettyPic,'D:\Imaging\2Color\PrettyPicAll','-dpdf');

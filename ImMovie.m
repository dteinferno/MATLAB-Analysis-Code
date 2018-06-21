function ImMovie(movFilename, tiffFilename, behavFilename, movMax, movSpeed)
% Create a movie that combines the view of the fly, the view of the Ca2+
% activity, the fly's heading, and the fly's position in the virtual world.
%
% ARGUMENTS:
%    INPUTS:
%    movFilename: .seq file containing the movie of the fly walking on the
%    ball
%    tiffFilename: .tif file containing the Ca2+ imaging stack
%    behavFilename: .txt file containing the behavioral information
%    movMax: multiplier for the maximum of the Ca2+ stack
%    movSpeed: multiplier for the playback speed

    % Load the movie
    [headerInfo, imgOut] = Norpix2MATLAB_All(movFilename);
    
    % Load the imaging stack and get the behavioral data
%     [fullpath, stackMaxInt, stackMean, positionDat] = VRImDatLoadDiscardNoPrompt(tiffFilename,strcat(pwd,'\'),behavFilename,strcat(pwd,'\'),1);
    [stackMaxInt, stackMean] = ImDatLoadBigtiff(tiffFilename,strcat(pwd,'\'));
    positionDat = VRDatLoad(behavFilename,strcat(pwd,'\'),0);
    
    % Gaussian filter the stack
    gaussianSize = [5 5];
    gaussianSigma = 2;
    Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
    stackXYfiltMax = double(zeros(size(stackMaxInt)));
    
    h = waitbar(0.0,'Gaussian filtering stack...');
    set(h,'Position',[50 50 360 72]);
    set(h,'Name','Gaussian filtering TIFF stack...');
    for i = 1:size(stackMaxInt,3)
        if mod(i,100)==0
            waitbar(i/length(stackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stackMaxInt))]);
        end
        stackXYfiltMax(:,:,i) = imfilter(stackMaxInt(:,:,i),Gxy,'replicate','conv');
    end
    delete(h);
    
    % Match the behavior to the imaging
    num_planes = round(length(positionDat.tFrameGrab)/size(stackMaxInt,3));
    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
    maxFG = min(round(length(positionDat.tFrameGrab)/num_planes),size(stackMaxInt,3));

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
                
    % Find the min and max of the imaging stack and the imaging rate
    minCa = min(min(min(stackXYfiltMax)));
    maxCa = movMax*max(max(max(stackXYfiltMax)));
    imFreq = 10000/mean(diff(positionDat.tFrameGrab))/num_planes;
    
    % Specify the length of past and future time to plot and create a
    % corresponding color scheme
    tPre = 1;
    tPre = round(tPre*imFreq);
    tPost = 1;
    tPost = round(tPost*imFreq);
    
    % general graphics, this will apply to any figure you open (groot is the default figure object).
    % I have this in my startup.m file, so I don't have to retype these things whenever plotting a new fig.
    set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 8, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(movFilename(1:end-4),'.avi'));
    writerObj.FrameRate= movSpeed*imFreq;
    open(writerObj);
    frameNow = figure('Position',[50 50 1000 1000]);
    
    % Generate the movie
    for fm = 1:maxFG-minFG+1
        subplot(2,2,1);
        cla;
        imagesc(flipud(fliplr(stackXYfiltMax(:,:,fm+minFG-1))),[minCa maxCa]);

        subplot(2,2,2);
        cla;
        hold on;
        if fm > tPre & fm < maxFG-minFG+1-tPost
            scatter(OffsetForMatch(fm-tPre:fm-1),OffsetLatMatch(fm-tPre:fm-1),10,[0.25 0.25 0.25],'filled');
            scatter(OffsetForMatch(fm+1:fm+tPost),OffsetLatMatch(fm+1:fm+tPost),10,[0.25 0.25 0.25]);
            scatter(OffsetForMatch(fm),OffsetLatMatch(fm),20,'k','filled');
            posNow = [OffsetForMatch(fm) OffsetLatMatch(fm)];
            line([posNow(1) posNow(1)+0.5*cos(OffsetRotMatch(fm,2))],...
                [posNow(2) posNow(2)+0.5*sin(OffsetRotMatch(fm,2))],...
                'Color','k');
        elseif fm > tPre+tPost
            scatter(OffsetForMatch(fm-tPre:fm-1),OffsetLatMatch(fm-tPre:fm-1),10,[0.25 0.25 0.25],'filled');
            scatter(OffsetForMatch(fm:end),OffsetLatMatch(fm:end),10,[0.25 0.25 0.25]);
            scatter(OffsetForMatch(fm),OffsetLatMatch(fm),20,'k','filled');
            posNow = [OffsetForMatch(fm) OffsetLatMatch(fm)];
            line([posNow(1) posNow(1)+0.5*cos(OffsetRotMatch(fm,2))],...
                [posNow(2) posNow(2)+0.5*sin(OffsetRotMatch(fm,2))],...
                'Color','k');
        else
            scatter(OffsetForMatch(1:fm),OffsetLatMatch(1:fm),10,[0.25 0.25 0.25],'filled');
            scatter(OffsetForMatch(fm+1:fm+tPost),OffsetLatMatch(fm+1:fm+tPost),10,[0.25 0.25 0.25]);
            scatter(OffsetForMatch(fm),OffsetLatMatch(fm),20,'k','filled');
            posNow = [OffsetForMatch(fm) OffsetLatMatch(fm)];
            line([posNow(1) posNow(1)+0.5*cos(OffsetRotMatch(fm,2))],...
                [posNow(2) posNow(2)+0.5*sin(OffsetRotMatch(fm,2))],...
                'Color','k');
        end
        axis equal;
        xlim([-10 10]);
        ylim([-10 10]);
        xlabel('x position (cm)');
        ylabel('y position (cm)');

        subplot(2,2,3);
        cla;
        imshow(imgOut(:,:,round((fm+minFG-1)*num_planes/2)),[0 255]);

        subplot(2,2,4);
        cla;
        hold on;
        if fm > tPre & fm < maxFG-minFG+1-tPost
            plot(OffsetRotMatch(fm-tPre:fm,1),OffsetRotMatch(fm-tPre:fm,2),'Color',[0.25 0.25 0.25]);
            plot(OffsetRotMatch(fm:fm+tPost,1),OffsetRotMatch(fm:fm+tPost,2),'Color',[0.25 0.25 0.25],'LineStyle','--');
            scatter(OffsetRotMatch(fm,1),OffsetRotMatch(fm,2),20,'k','filled');
        elseif fm > tPre+tPost
            plot(OffsetRotMatch(fm-tPre:fm,1),OffsetRotMatch(fm-tPre:fm,2),'Color',[0.25 0.25 0.25]);
            plot(OffsetRotMatch(fm:end,1),OffsetRotMatch(fm:end,2),'Color',[0.25 0.25 0.25],'LineStyle','--');
            scatter(OffsetRotMatch(fm,1),OffsetRotMatch(fm,2),20,'k','filled');
        else
            plot(OffsetRotMatch(1:fm,1),OffsetRotMatch(1:fm,2),'Color',[0.25 0.25 0.25]);
            plot(OffsetRotMatch(fm:fm+tPost,1),OffsetRotMatch(fm:fm+tPost,2),'Color',[0.25 0.25 0.25],'LineStyle','--');
            scatter(OffsetRotMatch(fm,1),OffsetRotMatch(fm,2),20,'k','filled');
        end
        axis tight;
        ylim([-pi pi]);
        xlabel('time (s)');
        ylabel('heading (rad)');
        
        frame = getframe(frameNow);
        writeVideo(writerObj,frame);
    end
    
    delete(frameNow);
    close(writerObj);
end
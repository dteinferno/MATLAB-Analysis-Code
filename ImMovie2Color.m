function ImMovie2Color(movFilename, tiffFilename, behavFilename, rotAng, movMax, movSpeed)
% Create a two-color movie that combines the view of the fly,
% the view of the Ca2+ activity, the fly's heading, and the fly's position
% in the virtual world.
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
    [RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] = ImDatLoadBigtiff2Color(tiffFilename,strcat(pwd,'\'),0);
    positionDat = VRDatLoad(behavFilename,strcat(pwd,'\'),0);
    
    % Subtract the background
%     RMean = mean(RstackMaxInt,3);
%     GMean = mean(GstackMaxInt,3);
%     for frame = 1:size(RstackMaxInt,3)
%         RstackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
%         GstackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
%     end
    RstackSub = RstackMaxInt;
    GstackSub = GstackMaxInt;
    
    % Gaussian filter the stacks
    gaussianSize = [5 5];
    gaussianSigma = 2;
    Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
    RstackXYfiltMax = double(zeros(size(RstackSub)));
    GstackXYfiltMax = double(zeros(size(RstackSub)));
    
    h = waitbar(0.0,'Gaussian filtering stack...');
    set(h,'Position',[50 50 360 72]);
    set(h,'Name','Gaussian filtering TIFF stack...');
    for i = 1:size(RstackMaxInt,3)
        if mod(i,100)==0
            waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
        end
        RstackXYfiltMax(:,:,i) = imfilter(RstackSub(:,:,i),Gxy,'replicate','conv');
        GstackXYfiltMax(:,:,i) = imfilter(GstackSub(:,:,i),Gxy,'replicate','conv');
    end
    delete(h);
    
    % Match the behavior to the imaging
    num_planes = round(length(positionDat.tFrameGrab)/size(RstackMaxInt,3));
    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
    maxFG = min(round(length(positionDat.tFrameGrab)/num_planes),size(RstackMaxInt,3));

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
    RmaxCa = movMax*max(max(max(RstackXYfiltMax)));
    GmaxCa = movMax*max(max(max(GstackXYfiltMax)));
    RminCa = min(min(min(RstackXYfiltMax)));% RmaxCa/10;
    GminCa = min(min(min(GstackXYfiltMax)));% GmaxCa/10;
    imFreq = 10000/mean(diff(positionDat.tFrameGrab))/num_planes; 
    
    % Open the movie object and make the figure
    writerObj = VideoWriter(strcat(movFilename(1:end-4),'.avi'));
    writerObj.FrameRate= movSpeed*imFreq;
    open(writerObj);
    frameNow = figure('Position',[50 50 1000 1000],'Color','k');
    
    % Generate the movie
    for fm = 1:maxFG-minFG+1
        
        % Create the colored images from the two stacks
        curRFrame = (imrotate(flipud(fliplr(RstackXYfiltMax(:,:,fm+minFG-1))),rotAng,'crop')-...
            RminCa)./(RmaxCa-RminCa);
        curGFrame = (imrotate(flipud(fliplr(GstackXYfiltMax(:,:,fm+minFG-1))),rotAng,'crop')-...
            GminCa)./(GmaxCa-GminCa);

        RIm = zeros([size(curRFrame) 3]);
        RIm(:,:,1) = curRFrame;
        RIm(:,:,3) = curRFrame;
        GIm = zeros([size(curGFrame) 3]);
        GIm(:,:,2) = curGFrame;
        
        subplot(2,2,1);
        cla;
        imshow(RIm);
        axis off;

        subplot(2,2,2);
        cla;
        imshow(GIm);
        axis off;

        
        % Show the fly's behavior
        curBehavFrame = imgOut(:,:,round((fm+minFG-1)*num_planes/4));
        
        subplot(2,2,3);
        cla;
        imshow(curBehavFrame(100:300,50:400) ,[0 255]);

        subplot(2,2,4,'color','k');
        cla;
        hold on;
        xPos = sin(OffsetRotMatch(fm,2));
        yPos = cos(OffsetRotMatch(fm,2));
        line([0 xPos], [0 yPos],'Color','w');
        rectangle('Position',[-1 -1 2 2],'Curvature',1,'EdgeColor','w');
        xlim([-2 2]);
        ylim([-2 2]);
        title('heading (rad)','Color','w');
        axis off;
        
        frame = getframe(frameNow);
        writeVideo(writerObj,frame);
    end
    
    delete(frameNow);
    close(writerObj);
end
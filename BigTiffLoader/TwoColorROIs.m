function [positionG, positionR] = TwoColorROIs(GstackMaxIntRot, RstackMaxIntRot,imRegion)

    % Initialize arrays to hold the ROI coordinates
    positionG = {};
    positionR = {};
    
    % Make some pretty colormaps
    colorBins = 64;
    Reds = 1-0.5*brewermap(colorBins, 'Greens')-0.5*brewermap(colorBins, 'Blues');
    Greens = 1-0.5*brewermap(colorBins, 'Reds')-0.5*brewermap(colorBins, 'Blues');

    % Normalize the mean images for plotting
    meanStack_G = mean(GstackMaxIntRot,3);
    meanStackColored_G = zeros([size(meanStack_G) ,3]);
    meanStackNorm_G = (meanStack_G-min(min(meanStack_G)))./(max(max(meanStack_G))-min(min(meanStack_G)));

    meanStack_R = mean(RstackMaxIntRot,3);
    meanStackColored_R = zeros([size(meanStack_R) ,3]);
    meanStackNorm_R = (meanStack_R-0.5*mean(mean(meanStack_R)))./(max(max(meanStack_R))-0.5*mean(mean(meanStack_R)));

    % Color the images with the colormaps
    for x = 1:size(meanStackNorm_G,1)
        for y = 1:size(meanStackNorm_G,2)
            mapBin_G = max(1,ceil(colorBins*meanStackNorm_G(x,y)));
            meanStackColored_G(x,y,:) = Greens(mapBin_G,:);

            mapBin_R = max(1,ceil(colorBins*meanStackNorm_R(x,y)));
            meanStackColored_R(x,y,:) = Reds(mapBin_R,:);
        end
    end
    
    % Plot the two colors
    hf = figure('units','normalized','outerposition',[0 0 1 1],'Color','k');
    subplot(2,2,1);
    image(meanStackColored_G);
    axis equal;
    axis off;

    subplot(2,2,3);
    image(meanStackColored_R);
    axis equal;
    axis off;

    subplot(2,2,[2 4]);
    image(meanStackColored_G + meanStackColored_R);
    axis equal;
    axis off;
    
    % Ask if the user wants to draw the same ref ROI for both images or
    % wants a separate one for each
    ss = input('Same (0) or separate (1) ROIs?');
    
    if ss == 0
        subplot(2,2,[2 4]);
        title('Draw ROIs here','Color','w');
        position = OneColorROIs(GstackMaxIntRot, imRegion, {})
        positionG = position;
        positionR = position;
        
        done = input('Fix the ROIs (y or n)');
        
        while strcmp(done,'y')
            subplot(2,2,1);
            cla;
            image(meanStackColored_G);
            axis equal;
            axis off;
            for pgs = 1:length(position)
                impoly(gca,position{pgs});
            end

            subplot(2,2,3);
            cla;
            image(meanStackColored_R);
            axis equal;
            axis off;
            for pgs = 1:length(position)
                impoly(gca,position{pgs});
            end
            
            subplot(2,2,[2 4]);
            cla;
            image(meanStackColored_G + meanStackColored_R);
            axis equal;
            axis off;
            
            position = OneColorROIs(GstackMaxIntRot, imRegion, position);
            positionG = position;
            positionR = position;
            
            done = input('Fix the ROIs (y or n)');
        end
    else
        for channelNow = 1:2
            if channelNow == 1
                stackMaxIntRot = GstackMaxIntRot;
                meanStackColored = meanStackColored_G;
            else
                stackMaxIntRot = RstackMaxIntRot;
                meanStackColored = meanStackColored_R;
            end
            
            subplot(2,2,2*channelNow-1);
            title('Draw ROIs here','Color','w');
            position = OneColorROIs(stackMaxIntRot, imRegion, {})
            
            done = input('Fix the ROIs (y or n)');

            while strcmp(done,'y')
                subplot(2,2,2*channelNow-1);
                cla;
                image(meanStackColored);
                axis equal;
                axis off;
                position = OneColorROIs(stackMaxIntRot, imRegion, position);
                
                done = input('Fix the ROIs (y or n)');
            end
                
            if channelNow == 1
                positionG = position;
            else
                positionR = position;
            end
        end
    end
    delete(hf);
end
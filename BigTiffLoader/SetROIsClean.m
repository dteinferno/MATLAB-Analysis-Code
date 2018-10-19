%% Get an ROI from just the first few volumes of a tif
function [ROIs] = SetROIsClean(allPathname,tifName,imRegion,numColors, stackMaxInt)
    
     % Rotate the stack if the EB was imaged
    if strcmp(imRegion{1},'EB') || strcmp(imRegion{1},'FB')
        rotAng = -70;
    else
        rotAng = 90;
    end
    if numColors == 1
        % Rotate the stacks so that the EB is aligned
        stackMean = mean(stackMaxInt{1},3);
        stackMean = imrotate(stackMean,rotAng);
        stackMaxIntRot = zeros([size(stackMean) size(stackMaxInt{1},3)]);
        for frame = 1:size(stackMaxIntRot,3)
            stackMaxIntRot(:,:,frame) = imrotate(stackMaxInt{1}(:,:,frame),rotAng);
        end
    elseif numColors == 2
        RMean = mean(stackMaxInt{2},3);
        RMean = imrotate(RMean,rotAng);
        GMean = mean(stackMaxInt{1},3);
        GMean = imrotate(GMean,rotAng);
        RstackMaxIntRot = zeros([size(RMean) size(stackMaxInt{2},3)]);
        GstackMaxIntRot = zeros([size(GMean) size(stackMaxInt{1},3)]);
        for frame = 1:size(RstackMaxIntRot,3)
            RstackMaxIntRot(:,:,frame) = imrotate(stackMaxInt{2}(:,:,frame),rotAng);
            GstackMaxIntRot(:,:,frame) = imrotate(stackMaxInt{1}(:,:,frame),rotAng);
        end
    end

    % Load the reference ROIs
    newPath = strrep(allPathname,'Downloads\Data','Documents\RawAnalysis');
    pathParts = strsplit(tifName,'_');
    IDLoc = strfind(tifName,pathParts{end-1});
    ROIfName = strcat(newPath,tifName(1:IDLoc-2));
    load(strcat(ROIfName,'_ReferenceROIs'));
    
    % Get the ROIs
    if numColors == 1
        stackSize = [size(stackMaxIntRot,1) size(stackMaxIntRot,2)];
    else
        stackSize = [size(RstackMaxIntRot,1) size(RstackMaxIntRot,2)];
    end
    
    if numColors == 1
        % Show the stack
        A = mean(stackMaxIntRot,3);
        hf = figure;
        hold;
        imshow(A,[0 max(max(A))]);
    
        ROIsNow = MakeROI(stackSize, imRegion{1}, position, A);
        ROIs{1} = ROIsNow;
        
        delete(hf);
    elseif numColors == 2
        
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
            ROIsNow = MakeROI(stackSize, imRegion{1}, positionG,meanStackColored_G + meanStackColored_R);
            ROIs{1} = ROIsNow;
            ROIs{2} = ROIsNow;
        else
            for channelNow = 1:2
                if channelNow == 1
                    position = positionG;
                    meanImage = meanStackColored_G;
                else
                    position = positionR;
                    meanImage = meanStackColored_R;
                end
                
                subplot(2,2,2*channelNow-1);
                title('Draw ROIs here','Color','w');
                ROIs = MakeROI(stackSize, imRegion{channelNow}, position, meanImage);
                
                if channelNow == 1
                    ROIs{1} = ROIsNow;
                else
                    ROIs{2} = ROIsNow;
                end
            end
        end
        delete(hf);
    end
end
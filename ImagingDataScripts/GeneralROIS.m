%% Function to calculate the intensity over time for arbitrary ROIs

function [ROIs, ROIaveREF] = GeneralROIS(stack)

    % Show the mean intensity and place some ROIs
    A = mean(stack,3);
    hf = figure('units','normalized','outerposition',[0 0 1 1]);
    hold;
    imshow(A,[0 max(max(A))]);

    num_ellipses = input('Number of ellipses?');
    num_pgons = input('Number of polygons?');
    num_ROIs = num_ellipses+num_pgons;

    ROIs = zeros(num_ROIs,size(stack,1),size(stack,2));
    for els = 1:num_ellipses
        h = imellipse;
        position = wait(h);
        ROIs(els,:,:) = roipoly(A,position(:,1),position(:,2));
    end
    for pgs = 1:num_pgons
        h = impoly;
        position = wait(h);
        ROIs(num_ellipses+pgs,:,:) = roipoly(A,position(:,1),position(:,2));
    end

    % Plot the ROI on each figure and calculate the average
    numFrames = size(stack,3);
    ROIaveREF = zeros(num_ROIs,numFrames);
    h = waitbar(0.0,'Calculating ROIs...');
    set(h,'Position',[50 50 360 72]);
    set(h,'Name','Calculating ROIs...');
    for i = 1:numFrames
        if mod(i,100)==0
            waitbar(i/numFrames,h,['Calculating frame# ' num2str(i) ' out of ' num2str(numFrames)]);
        end
        for incROI = 1:num_ROIs
            A = squeeze(stack(:,:,i));
            ROI = A(logical(squeeze(ROIs(incROI,:,:))));
            ROIaveREF(incROI,i) = mean2(ROI);
        end
    end
    delete(h);
    
    delete(hf);

end

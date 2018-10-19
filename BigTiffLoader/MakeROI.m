function ROIs = MakeROI(stackSize, imRegion, position, meanImage)
    
    % Place and make the ROIs
    if strcmp(imRegion,'EB')
        
        % Set the number of wedges
        numROIs = 16;

        % Initialize the ROI array
        ROIs = zeros([numROIs, stackSize]);
        
        % Find the change in angles
        deltaAng = 2*pi/numROIs;

        % Specify an outer ellipse
        hOut = impoly(gca,position{1});
        positionOut = wait(hOut);
        ellipseCenterOut = [mean(positionOut(:,1)) mean(positionOut(:,2))];
        ellipseHeightOut = max(positionOut(:,2)) - min(positionOut(:,2));
        ellipseWidthOut = max(positionOut(:,1)) - min(positionOut(:,1));

        % Specify an inner ellipse
        hIn = impoly(gca,position{2});
        positionIn = wait(hIn);
        ellipseCenterIn = [mean(positionIn(:,1)) mean(positionIn(:,2))];
        ellipseHeightIn = max(positionIn(:,2)) - min(positionIn(:,2));
        ellipseWidthIn = max(positionIn(:,1)) - min(positionIn(:,1));

        for incROI = 1:numROIs
            for angs=1:5
                xROI(angs) = ellipseCenterOut(1)+ellipseWidthOut/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
                yROI(angs) = ellipseCenterOut(2)+ellipseHeightOut/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-1));
            end
            for angs=6:10
                xROI(16-angs) = ellipseCenterIn(1)+ellipseWidthIn/2*cos(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
                yROI(16-angs) = ellipseCenterIn(2)+ellipseHeightIn/2*sin(deltaAng*(incROI-1)+deltaAng/4*(angs-6));
            end
            ROIs(incROI,:,:) = roipoly(meanImage, xROI,yROI);
        end
    else
        % Get the number of ROIs
        numROIs = length(position);

        % Initialize the ROI array
        ROIs = zeros([numROIs, stackSize]);

        % Place the ROIs
        for pgs = 1:numROIs
            pgons{pgs} = impoly(gca,position{pgs});
        end
        wait(pgons{1});

        % Get the ROIs
        for pgs = 1:numROIs
            posNow = getPosition(pgons{pgs});
            ROIs(pgs,:,:) = roipoly(meanImage,posNow(:,1),posNow(:,2));
        end
    end
end
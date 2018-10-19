function position = OneColorROIs(stackMaxIntRot, imRegion, position)

    % Get the ROI points
    if isempty(position)
        
        % For the EB, draw an inner and outer ellipse
        if strcmp(imRegion{1},'EB') 
            % Specify an outer ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            position{1} = posNow;

            % Specify an inner ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            position{2} = posNow;

        % For the PB, draw 18 glomeruli
        elseif strcmp(imRegion{1},'PB')

            % Initiate glomeruli in a PB like layout
            num_glomeruli = 18;
            pgons = {};
            scaleFact = round(size(stackMaxIntRot,1)/70);
            defaultGlomShape = [1,0;3,0;4,1;4,4;3,5;1,5;0,4;0,1];
            defaultGlomShape = scaleFact*defaultGlomShape;
            glomOffsets = [-7.5,2;-7.5,0.5;-6.5,-0.5;-5.5,-1.5;-4.5,-1.5;-3.5,-1.25;-2.5,-1;-1.5,-0.75;-0.5,-0.5];
            glomOffsets = vertcat(glomOffsets,flipud(glomOffsets.*[-1,1]));
            glomOffsets = 0.5*[size(stackMaxIntRot,1),size(stackMaxIntRot,2)]+...
                4*scaleFact*glomOffsets;

            % Give each of the glomeruli a different color to make them easier
            % to track
            for pgs = 1:num_glomeruli
                pgons{pgs} = impoly(gca,defaultGlomShape+glomOffsets(pgs,:));
                api = iptgetapi(pgons{pgs});
                api.setColor([(pgs-1)/(num_glomeruli-1) 1 1-(pgs-1)/(num_glomeruli-1)]);
            end

            wait(pgons{1});

            % Get the position of each glomerulus
            for pgs = 1:num_glomeruli
                posNow = getPosition(pgons{pgs});
                if contains(tifName,'Anatomy')
                    posNow = 256./size(stackMaxIntRot,1).*posNow;
                end

                position{pgs} = posNow;
            end

        % For the FB, draw 8 polygons around the 8 regions
        elseif strcmp(imRegion{1},'FB')
            for pgs = 1:8
                h = impoly;
                posNow = wait(h);
                position{pgs} = posNow;
            end

        % For the NO, draw an ellipse around each side
        elseif strcmp(imRegion{1},'NO')
            % Specify a left ellipse
            hOut = imellipse;
            posNow = wait(hOut);
            position{1} = posNow;

            % Specify a right ellipse
            hIn = imellipse;
            posNow = wait(hIn);
            position{2} = posNow;

        % For any other shapes, specify the number of ellipses and polygons to
        % draw
        elseif strcmp(imRegion{1},'other')
            num_ellipses = input('Number of ellipses?');
            num_pgons = input('Number of polygons?');
            num_ROIs = num_ellipses+num_pgons;

            ROIs = zeros(num_ROIs,size(stackMaxInt,1),size(stackMaxInt,2));
            ROIverts = cell(num_ROIs,1);
            for els = 1:num_ellipses
                h = imellipse;
                posNow = wait(h);
                position{els} = posNow;
            end
            for pgs = 1:num_pgons
                h = impoly;
                posNow = wait(h);
                position{num_ellipses+pgs} = posNow;
            end
        else
            disp('Functionality not yet enabled');
        end
        
    % Fix the ROI points
    else
        pgons = {};
        for pgs = 1:length(position)
            pgons{pgs} = impoly(gca,position{pgs});
        end
        wait(pgons{1});

        for pgs = 1:length(position)
            position{pgs} = getPosition(pgons{pgs});
        end
    
    end
end
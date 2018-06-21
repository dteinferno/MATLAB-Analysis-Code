%% Function to calculate the intensity over time for arbitrary ROIs

function [ROIs] = GeneralROIsNew(allPathname,tifName)

    % Specify the number of volumes to look at
    numVolumes = 100;

    % Get the 2P info;
    fullpath = strcat(allPathname,tifName);
    info = imfinfo(fullpath);
    hTiff = Tiff(fullpath);
    configStr = hTiff.getTag('Software');
    planeLoc = strfind(configStr, 'numFramesPerVolume');
    discardLoc = strfind(configStr, 'numDiscardFlybackFrames');
    num_planes = str2num(configStr(planeLoc+21:planeLoc+22));
    num_discards = str2num(configStr(discardLoc+26));
    width = info(1).Width;
    height = info(1).Height;

    % Find the max intensity projection
    stackMaxInt = double(zeros(height,width,numVolumes));
    miniStack = double(zeros(height,width,num_planes-num_discards));
    for incIm = 1:numVolumes*num_planes
        if mod(incIm,num_planes) ~= 0
            miniStack(:,:,1+mod(incIm-1,num_planes-num_discards)) = double(imread(fullpath, incIm, 'Info', info));
            if (1+mod(incIm-1,num_planes-num_discards) == num_planes-num_discards)
                stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
            end
        end
    end
    
    % Show the mean intensity and place some ROIs
    A = mean(stackMaxInt,3);
    hf = figure('units','normalized','outerposition',[0 0 1 1]);
    hold;
    imshow(A,[0 max(max(A))]);

    num_ellipses = input('Number of ellipses?');
    num_pgons = input('Number of polygons?');
    num_ROIs = num_ellipses+num_pgons;

    ROIs = zeros(num_ROIs,size(stackMaxInt,1),size(stackMaxInt,2));
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
    
    delete(hf);

end

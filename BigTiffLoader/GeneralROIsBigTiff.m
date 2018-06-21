%% Function to calculate the intensity over time for arbitrary ROIs

function [ROIs,ROIverts] = GeneralROIsBigTiff(allPathname,tifName,useROIs,prevROIs)
    % Create ROIs for processing a big tiff stack
    %
    % ARGUMENTS:
    %    INPUTS:
    %    allPathname: the directory where the .tif is located
    %    tiffFilename: .tif file containing the Ca2+ imaging stack
    %    useROIs: place previously generated ROIs on the mean stack
    %    prevROIs: filename of the previously generated ROIs

    % Default values for using previous ROIs
    if nargin < 3
        useROIs = 0;
        prevROIs = [];
    end
    
    % Get the tiff header info
    fullpath = strcat(allPathname,tifName);
    reader=ScanImageTiffReader(fullpath);
    desc = reader.metadata;
    planeLoc = strfind(desc, 'numFramesPerVolume');
    discardLoc = strfind(desc, 'numDiscardFlybackFrames');
    num_planes = str2num(desc(planeLoc+21:planeLoc+22));
    num_discards = str2num(desc(discardLoc+26));
    
    % Load the tifStack
    vol=ScanImageTiffReader(fullpath).data();
    width = size(vol,1);
    height = size(vol,2);
    numVolumes = floor(size(vol,3)/num_planes);

    % Find the max intensity projection
    stackMaxInt = double(zeros(height,width,numVolumes));
    miniStack = double(zeros(height,width,num_planes-num_discards));
    for incIm = 1:numVolumes*num_planes
        if mod(incIm,num_planes) ~= 0
            miniStack(:,:,1+mod(incIm-1,num_planes-num_discards)) = vol(:,:,incIm);
            if (1+mod(incIm-1,num_planes-num_discards) == num_planes-num_discards)
                stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
            end
        end
    end
    
    % Show the mean intensity
    A = mean(stackMaxInt,3);
    hf = figure('units','normalized','outerposition',[0 0 1 1]);
    hold;
    imshow(A',[0 max(max(A))]);
    ROIs = [];
    ROIverts = {};
    
    % Place new ROIs
    if ~useROIs
        num_ellipses = input('Number of ellipses?');
        num_pgons = input('Number of polygons?');
        num_ROIs = num_ellipses+num_pgons;

        ROIs = zeros(num_ROIs,size(stackMaxInt,1),size(stackMaxInt,2));
        ROIverts = cell(num_ROIs,1);
        for els = 1:num_ellipses
            h = imellipse;
            position = wait(h);
            ROIs(els,:,:) = roipoly(A,position(:,1),position(:,2));
            ROIverts{els}.X = position(:,1);
            ROIverts{els}.Y = position(:,2);
        end
        for pgs = 1:num_pgons
            h = impoly;
            position = wait(h);
            ROIs(num_ellipses+pgs,:,:) = roipoly(A,position(:,1),position(:,2));
            ROIverts{num_ellipses+pgs}.X = position(:,1);
            ROIverts{num_ellipses+pgs}.Y = position(:,2);
        end
    else % use previous ROIs
        num_ROIs = size(prevROIs,1);
        ROIs = zeros(num_ROIs,size(stackMaxInt,1),size(stackMaxInt,2));
        ROIverts = cell(num_ROIs,1);
        for plROI = 1:num_ROIs
            patch('XData',prevROIs{plROI}.X,'YData',prevROIs{plROI}.Y,'FaceColor','none','EdgeColor','w','LineWidth',0.1);
        end
        for plROI = 1:num_ROIs
            h = impoly(gca,horzcat(prevROIs{plROI}.X,prevROIs{plROI}.Y));
            position = wait(h);
            ROIs(plROI,:,:) = roipoly(A,position(:,1),position(:,2));
            ROIverts{plROI}.X = position(:,1);
            ROIverts{plROI}.Y = position(:,2);
        end
    end
    delete(hf);
    clear vol;
end
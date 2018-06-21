function ROIInfo = ROIExtract(stack)

% Show the mean intensity and place some ROIs
A = mean(stack,3);
hf = figure;
hold;
imshow(A,[0 max(max(A))]);


num_ellipses = input('Number of ellipses?');
num_pgons = input('Number of polygons?')
num_ROIs = num_ellipses+num_pgons;

ROIs = zeros(num_ROIs,size(stack,1),size(stack,2));
for els = 1:num_ellipses
    h = imellipse;
    position = wait(h);
    ROIs(els,:,:) = roipoly(A,position(:,1),position(:,2));
    ROIInfo(els).class = 'imellipse';
    ROIInfo(els).position = position;
    ROIInfo(els).logicalROI(:,:) = roipoly(A,position(:,1),position(:,2));
end
for pgs = 1:num_pgons
    h = impoly;
    position = wait(h);
    ROIs(num_ellipses+pgs,:,:) = roipoly(A,position(:,1),position(:,2));
    ROIInfo(pgs+num_ellipses).class = 'impoly';
    ROIInfo(pgs+num_ellipses).position = position;
    ROIInfo(pgs+num_ellipses).logicalROI(:,:) = roipoly(A,position(:,1),position(:,2));
end
    
% Plot the ROI on each figure and calculate the average
ROIaveREF = zeros(num_ROIs,length(stack));
h = waitbar(0.0,'Calculating ROIs...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Calculating ROIs...');
for i = 1:length(stack)
    if mod(i,100)==0
        waitbar(i/length(stack),h,['Calculating frame# ' num2str(i) ' out of ' num2str(length(stack))]);
    end
    for incROI = 1:num_ROIs
        A = squeeze(stack(:,:,i));
        ROI = A(logical(squeeze(ROIs(incROI,:,:))));
        ROIaveREF(incROI,i) = mean2(ROI);
    end
end
delete(h);

ROIave = zeros(num_ROIs,length(stack));
for incROI = 1:num_ROIs
    ROILow = sort(squeeze(ROIaveREF(incROI,:)));
    ROIave(incROI,:) = ROIaveREF(incROI,:)./sum(ROILow(1:floor(end/10)));
    ROIInfo(incROI).ROIave = squeeze(ROIave(incROI,:));
end

end
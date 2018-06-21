%% Subtract the background from a stack
function stackBGsub = BGSub(stack)

% Background subtraction
A = mean(stack,3);
h = figure;
ROIREF = roipoly(A/max(max(A)));
delete(h);

stackBGsub = double(zeros(size(stack)));

h = waitbar(0.0,'Background subtract...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Background subtract...');
for i = 1:length(stack)
    if mod(i,100)==0
        waitbar(i/length(stack),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stack))]);
    end
    A = stack(:,:,i);
    ROI = A(logical(ROIREF));
    stackBGsub(:,:,i) = stack(:,:,i) - mean2(ROI);
end
delete(h);
end

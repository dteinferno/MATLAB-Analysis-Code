%% Gaussian filter the stack
function stackXYfilt = GaussFilt(stack)

% Gaussian Filter and Background Subtraction
% Gaussian Filter
gaussianSize = [2 2];
gaussianSigma = 0.5;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
stackXYfilt = double(zeros(size(stack)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:length(stack)
    if mod(i,100)==0
        waitbar(i/length(stack),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stack))]);
    end
    stackXYfilt(:,:,i) = imfilter(stack(:,:,i),Gxy,'replicate','conv');
end
delete(h);

end
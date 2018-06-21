function stackReg = imReg(stack, xPx, yPx)

% Create a blank stack for the registered images
stackReg = zeros(size(stack));

% Get the region of interest for registering the image
figure;
imshow(mean(stack,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);
regShow = imrect;
regPos = wait(regShow);
xRng = round(regPos(1)-xPx):round(regPos(1)+regPos(3)+xPx);
yRng = round(regPos(2)-yPx):round(regPos(2)+regPos(4)+yPx);

% Specify the registration image and the area to be registered
ref_im = stack(yRng,xRng,round(size(stack,3)/2));
orig_res = [0 0 size(stack,1) size(stack,2)];
regArea = [xRng(1)+xPx yRng(1)+yPx length(xRng)-xPx*2+1 length(yRng)-yPx*2+1];

% relative offset of position of subimages
rect_offset = [(regArea(1)-xRng(1))
               (regArea(2)-yRng(1))];
           
h = waitbar(0.0,'Registering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');

% Register each image
for imNow = 1:size(stack,3)
    
    waitbar(imNow/size(stack,3),h,['Loading frame# ' num2str(imNow) ' out of ' num2str(size(stack,3))]);
        
    % Choose subregion of image
    current_im = stack(:,:,imNow);
    sub_EB = imcrop(current_im,regArea);

    % Do normalized cross correlation
    c = normxcorr2(sub_EB,ref_im);

    % offset found by correlation
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(sub_EB,2))
                   (ypeak-size(sub_EB,1))];

    % total offset
    offset = corr_offset;
    xoffset = offset(1);
    yoffset = offset(2);

    % coordinates for placing it in the registered stack
    xbegin = round(xoffset+1);
    xend   = round(xoffset+ size(stack,2)-2*xPx);
    ybegin = round(yoffset+1);
    yend   = round(yoffset+size(stack,1)-2*yPx);

    if ybegin <1 || yend > size(stack,1) || xbegin <1 || xend > size(stack,2)
        xbegin = xPx+1;
        ybegin = yPx+1;
        xend = (size(stack,2)-xPx);
        yend = (size(stack,1)-yPx);
    end
        
    % put it in the stack
    stackReg(ybegin:yend,xbegin:xend,imNow) = stack(xPx+1:(size(stack,1)-xPx),yPx+1:(size(stack,2)-yPx),imNow);
end

delete(h);

figure;
subplot(2,1,1)
imshow(mean(stack,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);
subplot(2,1,2)
imshow(mean(stackReg,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);


end

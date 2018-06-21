
% Create a blank stack for the registered images
stackReg = zeros(size(stack));

% Specify the region of interest
xRng = 100-10:185+10;
yRng = 80-10:175+10;
figure;
subplot(2,1,1);
imshow(stack(:,:,round(size(stack,3)/2)),[min(min(stack(:,:,1))) 0.5*max(max(stack(:,:,1)))]);
subplot(2,1,2);
imshow(stack(yRng,xRng,round(size(stack,3)/2)),[min(min(stack(:,:,1))) 0.5*max(max(stack(:,:,1)))]);

% Specify the registration image and the area to be registered
ref_im = stack(yRng,xRng,round(size(stack,3)/2));
orig_res = [0 0 size(stack,1) size(stack,2)];
regArea = [xRng(1)+10 yRng(1)+10 length(xRng)-19 length(yRng)-19];

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
    xend   = round(xoffset+ size(stack,2)-20);
    ybegin = round(yoffset+1);
    yend   = round(yoffset+size(stack,1)-20);

    if ybegin <1 || yend > size(stack,1) || xbegin <1 || xend > size(stack,2)
        xbegin = xPx+1;
        ybegin = yPx+1;
        xend = yPx+1;
        yend = square_EB(4)+square_EB(2);
    end
        
    % put it in the stack
    stackReg(ybegin:yend,xbegin:xend,imNow) = stack(11:(size(stack,1)-10),11:(size(stack,2)-10),imNow);
end


delete(h);

figure;
subplot(2,1,1)
imshow(mean(stack,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);
subplot(2,1,2)
imshow(mean(stackReg,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);

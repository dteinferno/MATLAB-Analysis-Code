function [RegStackReg,OtherStackReg, RegCoords] = imRegSimple2Color(RegStack, OtherStack,Px)

% Create blank stacks for the registered images
imWidth = size(RegStack,1);
imHeight = size(RegStack,2);
RegStackReg = zeros(size(RegStack));
OtherStackReg = zeros(size(OtherStack));
RegCoords = zeros(4,size(RegStack,3));

% Specify the registration image and the area to be registered
ref_im = mean(RegStack,3);
orig_res = [0 0 size(RegStack,1) size(RegStack,2)];
regArea = [Px Px size(RegStack,1)-Px*2+1 size(RegStack,2)-Px*2+1];

% relative offset of position of subimages
rect_offset = [Px
               Px];
           
h = waitbar(0.0,'Registering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Registering images...');

% Register each image
for imNow = 1:size(RegStack,3)
    
    waitbar(imNow/size(RegStack,3),h,['Loading frame# ' num2str(imNow) ' out of ' num2str(size(RegStack,3))]);
        
    % Choose subregion of image
    current_im = RegStack(:,:,imNow);
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
    xoffset = max(0,offset(1));
    yoffset = max(0,offset(2));

    % coordinates for placing it in the registered stack
    xbegin = round(xoffset+1);
    xend   = round(xoffset+ size(RegStack,2)-2*Px);
    ybegin = round(yoffset+1);
    yend   = round(yoffset+size(RegStack,1)-2*Px);

%     if ybegin <1 || yend > size(RegStack,1) || xbegin <1 || xend > size(RegStack,2)
%         xbegin = Px+1;
%         ybegin = Px+1;
%         xend = (size(RegStack,2)-Px);
%         yend = (size(RegStack,1)-Px);
%     end
        
    % put it in the stack
    RegStackRegTMP = zeros(imWidth+2*Px,imHeight+2*Px);
    RegStackRegTMP(ybegin:2*Px+yend,xbegin:2*Px+xend) = RegStack(:,:,imNow);
    RegStackReg(:,:,imNow) = imcrop(RegStackRegTMP,[Px+1 Px+1 imWidth-1 imHeight-1]);
    OtherStackRegTMP = zeros(imWidth+2*Px,imHeight+2*Px);
    OtherStackRegTMP(ybegin:2*Px+yend,xbegin:2*Px+xend) = OtherStack(:,:,imNow);
    OtherStackReg(:,:,imNow) = imcrop(OtherStackRegTMP,[Px+1 Px+1 imWidth-1 imHeight-1]);
    
    RegCoords(1,imNow) = xbegin;
    RegCoords(2,imNow) = xend;
    RegCoords(3,imNow) = ybegin;
    RegCoords(4,imNow) = yend;
end

delete(h);

% figure;
% subplot(2,1,1)
% imshow(mean(stack,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);
% subplot(2,1,2)
% imshow(mean(stackReg,3),[min(min(mean(stack,3))) max(max(mean(stack,3)))]);

end

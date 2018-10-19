function stackMaxInt = MaxIntBigTiff(allPathname,tifName,numColors)
    % Get the tiff header info
    fullpath = strcat(allPathname,tifName);
    reader = ScanImageTiffReader(fullpath);
    desc = reader.metadata;
    planeLoc = strfind(desc, 'numFramesPerVolume');
    discardLoc = strfind(desc, 'numDiscardFlybackFrames');
    num_planes = str2num(desc(planeLoc+21:planeLoc+22));
    num_discards = str2num(desc(discardLoc+26));

    % Load the tifStack
    vol=ScanImageTiffReader(fullpath).data();
    width = size(vol,1);
    height = size(vol,2);
    num_images = size(vol,3);
    numVolumes = floor(num_images/num_planes);
    if numColors == 2
        numVolumes = numVolumes/2;
    end
    
    % Find the max intensity projection
    if numColors == 1
        stackMaxInt{1} = double(zeros(height,width,numVolumes));
        miniStack = double(zeros(height,width,(num_planes-num_discards)));
        h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
        set(h,'Position',[50 50 360 72]);
        set(h,'Name','Loading TIFF stack...');
        for incIm = 1:num_images
            if mod(incIm,100)==0
                waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
            end
            if mod(ceil(incIm/2),num_planes) ~= 0
                miniStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                if (mod(incIm,num_planes) == num_planes-num_discards)
                    stackMaxInt{1}(:,:,ceil(incIm/(num_planes))) = max(miniStack,[],3);
                end
            end
        end
        delete(h);
    elseif numColors == 2
        stackMaxInt{1} = double(zeros(height,width,numVolumes));
        stackMaxInt{2} = double(zeros(height,width,numVolumes));
        RminiStack = double(zeros(height,width,(num_planes-num_discards)));
        GminiStack = double(zeros(height,width,(num_planes-num_discards)));
        h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
        set(h,'Position',[50 50 360 72]);
        set(h,'Name','Loading TIFF stack...');
        for incIm = 1:num_images
            if mod(incIm,100)==0
                waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
            end
            if mod(ceil(incIm/2),num_planes) ~= 0
                if mod(incIm,2)
                    RminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                else
                    GminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
                end
                if (mod(incIm/2,num_planes) == num_planes-num_discards)
                    stackMaxInt{1}(:,:,ceil(incIm/(2*num_planes))) = max(GminiStack,[],3);
                    stackMaxInt{2}(:,:,ceil(incIm/(2*num_planes))) = max(RminiStack,[],3);
                end
            end
        end
        delete(h);
    end
end
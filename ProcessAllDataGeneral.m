% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

% general graphics, this will apply to any figure you open (groot is the default figure object).
% I have this in my startup.m file, so I don't have to retype these things whenever plotting a new fig.
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 8, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);

% Specify the directory
allPathname = uigetdir('D:/Imaging','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

% Step through the files, and get an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifName
        ROIs = GeneralROIsNew(allPathname,tifName);
        save(strcat(tifName(1:end-4),'_ROIs','.mat'), 'ROIs');
        clear ROIs
    end
end

% Use the ROIs to get the DF/F for each ROI for the data
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifChunks = 1;
        % Extract the position and DF/F values
        for moveID = 3:length(fileNames)
            moveName =fileNames(moveID).name;
            tifNameParts = strsplit(tifName,'_');
            moveNameParts = strsplit(moveName,'_');
            % Find the associated position data 
            if (strcmpi(moveName(end-3:end),'.txt') &...
                    strcmp(tifName(end-5:end-4),moveName(end-5:end-4)) &...
                    strcmp(tifName(1:4),moveName(1:4)) &...
                    strcmp(moveNameParts{end-1},tifNameParts{end-1}))
                [fullpath, stackMaxInt, stackMean, positionDat] = VRImDatLoadDiscardNoPrompt(tifName,allPathname,moveName,allPathname,tifChunks);
                load(strcat(tifName(1:end-4),'_ROIs','.mat'));
                num_ROIs=size(ROIs,1);
%                 num_planes = length(positionDat.tFrameGrab)/length(stackMaxInt);
                num_planes = 10;
                
                 % Correct for movement
                stackRegMaxInt = stackMaxInt; %imRegSimple(stackMaxInt, 7);

                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                stackXYfiltMax = double(zeros(size(stackMaxInt)));
                
                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(stackMaxInt,3)
                    if mod(i,100)==0
                        waitbar(i/length(stackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(stackMaxInt))]);
                    end
                    stackXYfiltMax(:,:,i) = imfilter(stackRegMaxInt(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                ROIaveREFMax = zeros(num_ROIs,size(stackXYfiltMax,3));
                h = waitbar(0.0,'Calculating ROIs...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Calculating ROIs...');
                for i = 1:size(stackXYfiltMax,3)
                    if mod(i,10)==0
                        waitbar(i/size(stackXYfiltMax,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(stackXYfiltMax,3))]);
                    end
                    for incROI = 1:num_ROIs
                        A = squeeze(stackXYfiltMax(:,:,i));
                        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
                        ROIaveREFMax(incROI,i) = mean2(ROIMax);
                    end
                end
                delete(h);
                
                ROIaveMax = zeros(num_ROIs,length(stackMaxInt));
                for incROI = 1:num_ROIs
                    ROILow = sort(squeeze(ROIaveREFMax(incROI,:)));
                    ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./mean(ROILow(1:floor(end/10)));
                end
      
                save(strcat(tifName(1:end-4),'.mat'), 'fullpath','ROIaveMax','positionDat');
                
                if 0
                % Do some more preprocessing
                % Savitzky-Golay Filter the RFs
                sgolayOrder = 3;
                sgolayWindow = 11;
                ROIaveMaxFilt = sgolayfilt(ROIaveMax,sgolayOrder,sgolayWindow,[],2);

                % Match the behavior to the imaging
                num_planes = round(length(positionDat.tFrameGrab)/length(ROIaveMax));
                tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
                minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
                maxFG = min(round(length(positionDat.tFrameGrab)/num_planes),length(ROIaveMax));

                % Match the position data to the framegrab times
                OffsetRotMatch = zeros(maxFG-minFG+1,2);
                OffsetForMatch = zeros(maxFG-minFG+1,1);
                OffsetLatMatch = zeros(maxFG-minFG+1,1);
                for interp = minFG:maxFG
                    tMatch = find(positionDat.t >=...
                        (positionDat.t(1) +...
                        (positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                        positionDat.tVR(1))/10000));
                    OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
                    OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
                    OffsetForMatch(interp-minFG+1) = positionDat.OffsetFor(tMatch(1));
                    OffsetLatMatch(interp-minFG+1) = positionDat.OffsetLat(tMatch(1));
                end
                % Unwrap the rotation
                OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
                vRot = diff(OffsetRotUnwrap)./mean(diff(OffsetRotMatch(:,1)));
                vRot = sgolayfilt(vRot,sgolayOrder,sgolayWindow);
                vF = sqrt(diff(OffsetForMatch).^2+diff(OffsetLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));
                vF = sgolayfilt(vF,sgolayOrder,sgolayWindow);
                vSS = diff(OffsetLatMatch)./mean(diff(OffsetRotMatch(:,1)));
                
                % Compare the activity of each ROI to the stripe position,
                % rotational velocity, and forward velocity
                for incROI = 1:num_ROIs
                    % Calculate the circular mean and plot it on top of the activity
                    act = figure('units','normalized','outerposition',[0 0 1 1]);


                    % Plot the activity vs. various parameters
                    subplot(3,1,1);
                    hold on;
                    plot(OffsetRotMatch(:,1),ROIaveMaxFilt(incROI,minFG:maxFG)-1,'g');
                    plot(OffsetRotMatch(:,1),OffsetRotMatch(:,2)./(2*pi),'b');
                    set(gca,'FontSize',16);
                    ylabel('DF/F vs. position');
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);

                    subplot(3,1,2);
                    hold on;
                    plot(OffsetRotMatch(:,1),ROIaveMaxFilt(incROI,minFG:maxFG)-1,'g');
                    plot(OffsetRotMatch(1:end-1,1),vRot./max(vRot),'r');
                    set(gca,'FontSize',16);
                    ylabel('DF/F vs. vRot');
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);


                    subplot(3,1,3);
                    hold on;
                    plot(OffsetRotMatch(:,1),ROIaveMaxFilt(incROI,minFG:maxFG)-1,'g');
                    plot(OffsetRotMatch(1:end-1,1),vF./max(vF),'r');
                    set(gca,'FontSize',16);
                    ylabel('DF/F vs. vF');
                    xlabel('time (s)');
                    xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);

                    set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                    print(act,strcat(allPathname,tifName(1:end-4),'_Act',num2str(incROI)),'-dpdf');

                    delete(act);
                end
                end
                clear fullpath stackMean stackMaxInt positionDat ROIaveMax;
            end
        end
    end
end    

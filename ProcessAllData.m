%% Specify the directory
allPathname = uigetdir('D:/Imaging','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

%% Step through the files, and get an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-3:end),'.tif') & isempty(strfind(tifName,'Anatomy'))
        tifName
        ROIs = GetROIs(allPathname,tifName);
        save(strcat(tifName(1:end-4),'_ROIs','.mat'), 'ROIs');
        clear ROIs
    end
end

%% Use the ROIs to get the DF/F for each ROI for the data
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
                    strcmp(moveNameParts{end-1},tifNameParts{end-1}) &...
                    strcmp(moveNameParts{end-2},tifNameParts{end-2}))
                moveName
                [stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathname,0);
                positionDat = VRDatLoad(moveName,allPathname,0);
                
                load(strcat(tifName(1:end-4),'_ROIs','.mat'));
                num_ROIs=size(ROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(stackMaxInt);
                
                 % Correct for movement
                stackRegMaxInt = imRegSimple(stackMaxInt, 10);

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
                
                ROIaveMax = zeros(num_ROIs,size(stackMaxInt,3));
                for incROI = 1:num_ROIs
                    ROILow = sort(squeeze(ROIaveREFMax(incROI,:)));
                    ROIaveMax(incROI,:) = ROIaveREFMax(incROI,:)./mean(ROILow(1:floor(end/10)));
                end
      
                save(strcat(tifName(1:end-4),'.mat'), 'ROIaveMax','positionDat');
                
                % Do some more preprocessing
                % Savitzky-Golay Filter the RFs
                sgolayOrder = 3;
                sgolayWindow = 11;
                ROIaveMaxFilt = sgolayfilt(ROIaveMax,sgolayOrder,sgolayWindow,[],2);

                % Match the behavior to the imaging
                num_planes = round(length(positionDat.tFrameGrab)/size(ROIaveMax,2));
                tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
                minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
                maxFG = round(length(positionDat.tFrameGrab)/num_planes);

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
                
                % Calculate the circular mean and plot it on top of the activity
                act = figure('units','normalized','outerposition',[0 0 1 1]);

                % Calculate and plot the circular mean
                angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
                angsraw = angsraw';
                clear meanAngRaw;
                clear meanIntRaw;
                for ts = 1:size(ROIaveMax,2)
                    meanAngRaw(ts) = circ_mean(angsraw, squeeze(ROIaveMaxFilt(:,ts)));
                    meanIntRaw(ts) = circ_r(angsraw, squeeze(ROIaveMaxFilt(:,ts)));
                end

                meanAngRawPlot = meanAngRaw;
                for i=1:size(meanAngRaw,3)-1
                    if abs(meanAngRaw(i)-meanAngRaw(i+1)) > pi
                        meanAngRawPlot(i) = NaN;
                    end
                end

                % Plot the activity and circular mean
                subplot(3,1,1);
                imagesc(flipud(ROIaveMaxFilt(:,minFG:maxFG)));
                colormap(brewermap(64, 'Blues'));
                caxis([min(min(ROIaveMaxFilt)) max(max(ROIaveMaxFilt))]);
                set(gca,'FontSize',16);
                ylabel('position (rad)');
                title('maximum intensity');
                
                subplot(3,1,2);
                plot(OffsetRotMatch(:,1),meanAngRawPlot(minFG:maxFG),'b');
                set(gca,'FontSize',16);
                ylabel('PVA position (rad)');
                xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);
                
                subplot(3,1,3);
                plot(OffsetRotMatch(:,1),OffsetRotMatch(:,2),'k');
                set(gca,'FontSize',16);
                ylabel('fly heading (rad)');
                xlabel('time (sec)');
                xlim([OffsetRotMatch(1,1) OffsetRotMatch(end,1)]);

                set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                print(act,strcat(allPathname,tifName(1:end-4),'_BumpAct'),'-dpdf');
                
                delete(act);

                clear stackMean stackMaxInt positionDat ROIaveMax;
            end
        end
    end
end    

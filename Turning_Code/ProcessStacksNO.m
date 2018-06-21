% Specify the directory
allPathname = uigetdir('D:/Imaging','Select the directory');
fileNames = dir(allPathname);
cd(allPathname);
allPathname = strcat(allPathname,'\');

% Step through the files, and get an ROI for each set of tifs
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-4:end),'1.tif')
        tifName
        ROIs = GetROIsOnePlane(allPathname,tifName,6);
        save(strcat(tifName(1:end-10),'_NOROIs','.mat'), 'ROIs');
        clear ROIs
    end
end

% Use the ROIs to get the DF/F for each ROI for the data
for tifID = 3:length(fileNames)
    tifName = fileNames(tifID).name;
    % Process each trial
    if strcmp(tifName(end-4:end),'1.tif')
        tifChunks = 1;
        % Find how many tifs are in each stack
        for tifID2 = 3:length(fileNames)
            tifName2 = fileNames(tifID2).name;
            if strcmp(tifName(1:end-5),tifName2(1:end-5))
                tifChunks = max(tifChunks,str2num(tifName2(end-4)));
            end
        end
        % Extract the position and DF/F values
        for moveID = 3:length(fileNames)
            moveName =fileNames(moveID).name;
            tifNameParts = strsplit(tifName,'_');
            moveNameParts = strsplit(moveName,'_');
            % Find the associated position data 
            if (strcmpi(moveName(end-3:end),'.txt') &...
                    strcmp(tifName(end-11:end-10),moveName(end-5:end-4)) &...
                    strcmp(tifName(1:4),moveName(1:4)) &...
                    strcmp(moveNameParts{end-1},tifNameParts{end-2}))
                [fullpath, Rstack, Gstack, positionDat] = VRImDatLoadDiscardNoPromptOnePlane(tifName,allPathname,moveName,allPathname,tifChunks,6);
                load(strcat(tifName(1:end-10),'_NOROIs','.mat'));
                num_ROIs=size(ROIs,1);
                num_planes = length(positionDat.tFrameGrab)/length(Rstack);
                
                % Gaussian Filter
                gaussianSize = [2 2];
                gaussianSigma = 0.5;
                Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);

                RstackXYfilt = double(zeros(size(Rstack)));
                GstackXYfilt = double(zeros(size(Gstack)));
                
                
                h = waitbar(0.0,'Gaussian filtering stack...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Gaussian filtering TIFF stack...');
                for i = 1:size(Rstack,3)
                    if mod(i,100)==0
                        waitbar(i/length(Rstack),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(Rstack))]);
                    end
                    RstackXYfilt(:,:,i) = imfilter(Rstack(:,:,i),Gxy,'replicate','conv');
                    GstackXYfilt(:,:,i) = imfilter(Gstack(:,:,i),Gxy,'replicate','conv');
                end
                delete(h);
                
                % Plot the ROI on each figure and calculate the average for both the
                % Maximum Intesity and Mean projections
                RROIaveREF = zeros(num_ROIs,size(RstackXYfilt,3));
                GROIaveREF = zeros(num_ROIs,size(GstackXYfilt,3));
                h = waitbar(0.0,'Calculating ROIs...');
                set(h,'Position',[50 50 360 72]);
                set(h,'Name','Calculating ROIs...');
                for i = 1:size(RstackXYfilt,3)
                    if mod(i,10)==0
                        waitbar(i/size(RstackXYfilt,3),h,['Calculating frame# ' num2str(i) ' out of ' num2str(size(RstackXYfilt,3))]);
                    end
                    for incROI = 1:num_ROIs
                        A = squeeze(RstackXYfilt(:,:,i));
                        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
                        RROIaveREF(incROI,i) = mean2(ROIMax);
                        
                        A = squeeze(GstackXYfilt(:,:,i));
                        ROIMax = A(logical(squeeze(ROIs(incROI,:,:))));
                        GROIaveREF(incROI,i) = mean2(ROIMax);
                    end
                end
                delete(h);
                
                RROIave = zeros(num_ROIs,length(Rstack));
                GROIave = zeros(num_ROIs,length(Gstack));
                for incROI = 1:num_ROIs
                    ROILow = sort(squeeze(RROIaveREF(incROI,:)));
                    RROIave(incROI,:) = RROIaveREF(incROI,:)./mean(abs(ROILow(1:floor(end/10))));
                    
                    ROILow = sort(squeeze(GROIaveREF(incROI,:)));
                    GROIave(incROI,:) = GROIaveREF(incROI,:)./mean(abs(ROILow(1:floor(end/10))));
                end
      
                save(strcat(tifName(1:end-10),'_NO.mat'), 'fullpath','RROIave','GROIave','positionDat');
%                 save(strcat(tifName(1:end-10),'_NORaw.mat'), 'fullpath','RROIaveREF','GROIaveREF','positionDat');
                
%                 % Savitzky-Golay Filter the RFs
%                 sgolayOrder = 3;
%                 sgolayWindow = 11;
%                 RROIaveFilt = sgolayfilt(RROIave,sgolayOrder,sgolayWindow,[],2);
%                 GROIaveFilt = sgolayfilt(GROIave,sgolayOrder,sgolayWindow,[],2);
%                           
%                 % Match the behavior to the imaging
%                 tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
%                 minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
%                 maxFG = round(length(positionDat.tFrameGrab)/num_planes);
% 
%                 % Match the position data to the framegrab times
%                 OffsetRotMatch = zeros(maxFG-minFG+1,2);
%                 OffsetForMatch = zeros(maxFG-minFG+1,1);
%                 OffsetLatMatch = zeros(maxFG-minFG+1,1);
%                 for interp = minFG:maxFG
%                     tMatch = find(positionDat.t >=...
%                         (positionDat.t(1) +...
%                         (positionDat.tFrameGrab((interp-1)*num_planes+1)-...
%                         positionDat.tVR(1))/10000));
%                     OffsetRotMatch(interp-minFG+1,1) = positionDat.t(tMatch(1));
%                     OffsetRotMatch(interp-minFG+1,2) = pi/180*positionDat.OffsetRot(tMatch(1));
%                     OffsetForMatch(interp-minFG+1) = positionDat.OffsetFor(tMatch(1));
%                     OffsetLatMatch(interp-minFG+1) = positionDat.OffsetLat(tMatch(1));
%                 end
%                 % Unwrap the rotation
%                 OffsetRotUnwrap = UnWrap(OffsetRotMatch(:,2),2,0);
%                 vRot = diff(OffsetRotUnwrap)./mean(diff(OffsetRotMatch(:,1)));
%                 vRot = sgolayfilt(vRot,sgolayOrder,sgolayWindow);
%             
%                 % Plot the profiles
%                 act = figure('units','normalized','outerposition',[0 0 1 1]);
% 
%                 % Plot the activity
%                 subplot(4,1,1);
%                 plot(OffsetRotMatch(1:end-1,1),vRot,'k');
%                 set(gca,'FontSize',14);
%                 ylabel('v_{R} (rad/sec)');
%                 xlim([0 OffsetRotMatch(end,1)]);
%                 ylim([-2 2]);
%                 
%                 subplot(4,1,2);
%                 hold on;
%                 plot(OffsetRotMatch(:,1),RROIaveFilt(1,minFG:maxFG),'Color',[1 0.75 0]);
%                 plot(OffsetRotMatch(:,1),RROIaveFilt(2,minFG:maxFG),'Color',[1 0 0.75]);
%                 legend({'left','right'});
%                 set(gca,'FontSize',14);
%                 ylabel('\DeltaF/F');
%                 xlim([0 OffsetRotMatch(end,1)]);
%                 
%                 subplot(4,1,3);
%                 hold on;
%                 plot(OffsetRotMatch(:,1),GROIaveFilt(1,minFG:maxFG),'Color',[0.75 1 0]);
%                 plot(OffsetRotMatch(:,1),GROIaveFilt(2,minFG:maxFG),'Color',[0 1 0.75]);
%                 legend({'left','right'});
%                 set(gca,'FontSize',14);
%                 ylabel('\DeltaF/F');
%                 xlim([0 OffsetRotMatch(end,1)]);
%                 
%                 subplot(4,1,4);
%                 hold on;
%                 plot(OffsetRotMatch(1:end-1,1),vRot./max(abs(vRot)),'k');
%                 diffSig = GROIaveFilt(1,minFG:maxFG)-GROIaveFilt(2,minFG:maxFG);
%                 plot(OffsetRotMatch(:,1),diffSig./max(abs(diffSig)),'Color',[0 1 0]);
%                 set(gca,'FontSize',14);
%                 xlim([0 OffsetRotMatch(end,1)]);
%                 xlabel('time (sec)');
%                 
%                 set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%                 print(act,strcat(tifName(1:end-4),'_BumpActNOAlone'),'-dpdf'); 
%                 
%                 delete(act);

                clear fullpath Rstack Gstack positionDat RROIave GROIave;
                
            end
        end
    end
end    
%% Create array of bump shapes as sorted by velocity and sign of acceleration (for dark case)
% Time over which to smooth the turns
tSmooth = 1;

% Threshold velocities to consider
% velThreshLC = [0:20:60];
% velThreshUC = [20:20:80];
velThreshLC = [0:0.1:0.5];
velThreshUC = [0.1:0.1:0.6];
XCVals = {};
for flyID = 1:length(allFlyData)
    XCVals{flyID}.ID = allFlyData{flyID}.ID;
end

% Pull out the relevant cross sections for the dark cases.
for flyID = 1:length(allFlyData)
    
    % Pull out all trace for the dark trials
    for darkID = 1:length(allFlyData{flyID}.dark)
        if ~isempty(allFlyData{flyID}.dark{darkID})
            % Savitzky-Golay Filter the RFs
            sgolayOrder = 3;
            sgolayWindow = 7;
            RROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMax,sgolayOrder,sgolayWindow,[],2);
            GROIaveMaxFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMax,sgolayOrder,sgolayWindow,[],2);
            RROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.RROIaveMean,sgolayOrder,sgolayWindow,[],2);
            GROIaveMeanFilt = sgolayfilt(allFlyData{flyID}.dark{darkID}.GROIaveMean,sgolayOrder,sgolayWindow,[],2);

            % Match the behavior to the imaging
            num_planes = round(length(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab)/...
                length(allFlyData{flyID}.dark{darkID}.RROIaveMax));
            tSpan = allFlyData{flyID}.dark{darkID}.positionDat.t - ...
                allFlyData{flyID}.dark{darkID}.positionDat.t(1) + ...
                allFlyData{flyID}.dark{darkID}.positionDat.tVR(1)/10000;
            minFG = ceil(min(find(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab >=...
                allFlyData{flyID}.dark{darkID}.positionDat.tVR(1)))/num_planes);
            maxFG = round(length(allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab)/num_planes);

            % Match the position data to the framegrab times
            OffsetRotMatch = zeros(maxFG-minFG+1,2);
            OffsetForMatch = zeros(maxFG-minFG+1,1);
            OffsetLatMatch = zeros(maxFG-minFG+1,1);
            for interp = minFG:maxFG
                tMatch = find(allFlyData{flyID}.dark{darkID}.positionDat.t >=...
                    (allFlyData{flyID}.dark{darkID}.positionDat.t(1) +...
                    (allFlyData{flyID}.dark{darkID}.positionDat.tFrameGrab((interp-1)*num_planes+1)-...
                    allFlyData{flyID}.dark{darkID}.positionDat.tVR(1))/10000));
                OffsetRotMatch(interp-minFG+1,1) = allFlyData{flyID}.dark{darkID}.positionDat.t(tMatch(1));
                OffsetRotMatch(interp-minFG+1,2) = pi/180*allFlyData{flyID}.dark{darkID}.positionDat.OffsetRot(tMatch(1));
                OffsetForMatch(interp-minFG+1) = allFlyData{flyID}.dark{darkID}.positionDat.OffsetFor(tMatch(1));
                OffsetLatMatch(interp-minFG+1) = allFlyData{flyID}.dark{darkID}.positionDat.OffsetLat(tMatch(1));
            end
            
            OffsetForMatch = smooth(OffsetForMatch,round(tSmooth/mean(diff(OffsetRotMatch(:,1)))));
            
                    
            for velCut = 1:length(velThreshLC)
                % Sort into right and left turn bins by increasing or
                % decreasing velocity
                RturnsGinc = [];
                RturnsGdec = [];
                RturnsRinc = [];
                RturnsRdec = [];


                velThreshUCNow = velThreshUC(velCut).*mean(diff(OffsetRotMatch(:,1)));
                velThreshLCNow = velThreshLC(velCut).*mean(diff(OffsetRotMatch(:,1)));
                aFly = diff(diff(OffsetForMatch));

                for i=minFG+1:maxFG-2
                    if aFly(i+1-minFG) > 0 
                         if OffsetForMatch(i+1-minFG) - OffsetForMatch(i-minFG) > velThreshLCNow ...
                        & OffsetForMatch(i+1-minFG) - OffsetForMatch(i-minFG) < velThreshUCNow
                            RturnsGinc = horzcat(RturnsGinc, GROIaveMaxFilt(:,i));
                            RturnsRinc = horzcat(RturnsRinc, RROIaveMaxFilt(:,i));
                         end
                    elseif aFly(i+1-minFG) < 0
                        if OffsetForMatch(i+1-minFG) - OffsetForMatch(i-minFG) > velThreshLCNow ...
                                & OffsetForMatch(i+1-minFG) - OffsetForMatch(i-minFG) < velThreshUCNow
                            RturnsGdec = horzcat(RturnsGdec, GROIaveMaxFilt(:,i));
                            RturnsRdec = horzcat(RturnsRdec, RROIaveMaxFilt(:,i));
                        end
                    end
                end

                % Save the values to the array
                XCVals{flyID}.vel{velCut}.threshold = velThreshLC(velCut);
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsGinc = RturnsGinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsGdec = RturnsGdec;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsRinc = RturnsRinc;
                XCVals{flyID}.vel{velCut}.dark{darkID}.RturnsRdec = RturnsRdec;
            end
            
            clear RROIaveMaxFilt GROIaveMaxFilt OffsetRotMatch;
        end
    end
end

% Pool together the turns across trials
num_ROIs = 16;
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       RturnsGallGalignInc = [];
       RturnsGallGalignDec = [];
       RturnsRallGalignInc = [];
       RturnsRallGalignDec = [];
       
       RturnsGallRalignInc = [];
       RturnsGallRalignDec = [];
       RturnsRallRalignInc = [];
       RturnsRallRalignDec = [];
        for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
            RturnsGGalignInc = [];
            RturnsGGalignDec = [];
            RturnsRGalignInc = [];
            RturnsRGalignDec = [];
            
            RturnsGRalignInc = [];
            RturnsGRalignDec = [];
            RturnsRRalignInc = [];
            RturnsRRalignDec = [];
            
            % align the peaks by the green channel
            if ~isempty(XCVals{flyID}.vel{velCut}.dark{trial})
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i)));
                    RturnsGGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGinc(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignInc(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRinc(:,i),num_ROIs/2-pkPosRalign);
                end
                for i=1:size(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec,2)
                    pkPosRalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i)));
                    pkPosGalign = find(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i)==...
                        max(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i)));
                    RturnsGGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsRGalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i),num_ROIs/2-pkPosGalign);
                    RturnsGRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsGdec(:,i),num_ROIs/2-pkPosRalign);
                    RturnsRRalignDec(:,i) = circshift(XCVals{flyID}.vel{velCut}.dark{trial}.RturnsRdec(:,i),num_ROIs/2-pkPosRalign);
                end
            end
            
            RturnsGallGalignInc = horzcat(RturnsGallGalignInc,RturnsGGalignInc);
            RturnsRallGalignInc = horzcat(RturnsRallGalignInc,RturnsRGalignInc);
            RturnsGallGalignDec = horzcat(RturnsGallGalignDec,RturnsGGalignDec);
            RturnsRallGalignDec = horzcat(RturnsRallGalignDec,RturnsRGalignDec);
            
            RturnsGallRalignInc = horzcat(RturnsGallRalignInc,RturnsGRalignInc);
            RturnsRallRalignInc = horzcat(RturnsRallRalignInc,RturnsRRalignInc);
            RturnsGallRalignDec = horzcat(RturnsGallRalignDec,RturnsGRalignDec);
            RturnsRallRalignDec = horzcat(RturnsRallRalignDec,RturnsRRalignDec);

       end
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc = RturnsGallGalignInc;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc = RturnsRallGalignInc;
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec = RturnsGallGalignDec;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec = RturnsRallGalignDec;
       
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc = RturnsGallRalignInc;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc = RturnsRallRalignInc;
       XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec = RturnsGallRalignDec;
       XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec = RturnsRallRalignDec;
    end 
end

bumpCompRinc = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));
           
           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
           set(gca,'FontSize',14);
       end
       if flyID == 1
           title(strcat('v_{thresh} = ',num2str(XCVals{flyID}.vel{velCut}.threshold),'deg/sec'));
       end
       if flyID == numFlies
           xlabel('ROI');
       end
    end
    subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+1);
    ylabel(allFlyData{flyID}.ID);
end


bumpCompRdec = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID = 1:length(allFlyData)
    for velCut = 1:length(velThreshLC)       
       for trial=1:length(XCVals{flyID}.vel{velCut}.dark)
           RturnsG = XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec;
           RturnsR = XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec;
           RnetR = (RturnsR-min(median(RturnsR,2)))./(max(median(RturnsR,2))-min(median(RturnsR,2)));
           RnetG = (RturnsG-min(median(RturnsG,2)))./(max(median(RturnsG,2))-min(median(RturnsG,2)));

           subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+velCut);
           hold on;
           
           plot(quantile(RnetR,0.25,2),'b','LineWidth',1);
           plot(quantile(RnetR,0.75,2),'b','LineWidth',1);
           plot(median(RnetR,2),'b','LineWidth',4);
           text(15,1.1,num2str(size(RnetR,2)),'Color','b');
           plot(median(RnetG,2),'--g','LineWidth',1);
           
           xlim([1 16]);
           ylim([0 1.5]);
           set(gca,'FontSize',14);
       end
       if flyID == 1
           title(strcat('v_{thresh} = ',num2str(XCVals{flyID}.vel{velCut}.threshold),'deg/sec'));
       end
       if flyID == numFlies
           xlabel('ROI');
       end
    end
    subplot(numFlies,length(velThreshLC),(flyID-1)*length(velThreshLC)+1);
    ylabel(allFlyData{flyID}.ID);
end

% set(bumpCompRinc,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompRinc,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccInc'),'-dpdf');
% set(bumpCompRdec,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpCompRdec,strcat('Results\BumpCompGAlign_Red',redLine,'_Green',greenLine,'AccDec'),'-dpdf');

%% Get the difference between PVA direction for each fly
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=1:numFlies
    for velCut = 1:length(velThreshLC)       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec);
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec);
       
%        RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
%            XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
%            flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),...
%            flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec));
%        RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
%            XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
%            flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),...
%            flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec));

       % Calculate the PVA differences
       num_ROIs = 16;
       angsraw = (1:num_ROIs)*2*pi/num_ROIs+pi;
       angsraw = angsraw';
       PVAdiff(flyID,velCut) = circ_mean(-angsraw,median(RturnsR,2))+circ_mean(angsraw,median(RturnsG,2));
       PVAdiff(flyID,velCut) = min(PVAdiff(flyID,velCut),2*pi-PVAdiff(flyID,velCut));
       
    end

end

PVAdiff = PVAdiff;

figure;
hold on;
scatter(zeros(numFlies,1),180/pi*PVAdiff(:,1),'filled','k');
line([-0.1 0.1],[180/pi*mean(PVAdiff(:,1)) 180/pi*mean(PVAdiff(:,1))],'Color','k','LineWidth',2);
scatter(zeros(numFlies,1)+1,180/pi*PVAdiff(:,2),'filled','k');
line([0.9 1.1],[180/pi*mean(PVAdiff(:,2)) 180/pi*mean(PVAdiff(:,2))],'Color','k','LineWidth',2);
scatter(zeros(numFlies,1)+2,180/pi*PVAdiff(:,3),'filled','k');
line([1.9 2.1],[180/pi*mean(PVAdiff(:,3)) 180/pi*mean(PVAdiff(:,3))],'Color','k','LineWidth',2);
scatter(zeros(numFlies,1)+3,180/pi*PVAdiff(:,4),'filled','k');
line([2.9 3.1],[180/pi*mean(PVAdiff(:,4)) 180/pi*mean(PVAdiff(:,4))],'Color','k','LineWidth',2);
scatter(zeros(numFlies,1)+4,180/pi*PVAdiff(:,5),'filled','k');
line([3.9 4.1],[180/pi*mean(PVAdiff(:,5)) 180/pi*mean(PVAdiff(:,5))],'Color','k','LineWidth',2);
set(gca,'FontSize',14);
ylabel('PVA difference (deg)');
set(gca,'XTickLabel',{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'});
xlabel('vt bins (cm/sec)');
xlim([-0.5 4.5]);
% ylim([-10 40]);
% ylim([-5 30]);

% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,2),PVAdiff(:,3)),vertcat(zeros(10,1), ones(10,1), 2*ones(10,1)))
% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,2)),vertcat(zeros(5,1), ones(5,1)))
% anova1(vertcat(PVAdiff(:,2),PVAdiff(:,3)),vertcat(zeros(5,1), ones(5,1)))
% anova1(vertcat(PVAdiff(:,1),PVAdiff(:,3)),vertcat(zeros(5,1), ones(5,1)))

%% Get the peak heights and widths 
for flyID=[1 2 3 4 5 6 7 8 9 10]
    for velCut = 1:length(velThreshLC)       
       GturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec);
       GturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec);
       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec);
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec);

    if ~isempty(RturnsG) &  ~isempty(RturnsR) 
       % Calculate the PVA differences
       GpkMax(flyID, velCut) = max(median(RturnsG,2))-min(median(RturnsG,2));
       GpkRng = find(median(RturnsG,2)-min(median(RturnsG,2))>0.5*GpkMax(flyID, velCut));
       GpkWid(flyID, velCut) = GpkRng(end)-GpkRng(1);
    end
    if ~isempty(GturnsG) &  ~isempty(GturnsR)
       RpkMax(flyID, velCut) = max(median(GturnsR,2))-min(median(GturnsR,2));
       RpkRng = find(median(GturnsR,2)-min(median(GturnsR,2))>0.5*RpkMax(flyID, velCut));
       RpkWid(flyID, velCut) = RpkRng(end)-RpkRng(1);
       
    end
    end

end

figure;
subplot(2,2,1);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+0.1*i,squeeze(GpkMax(:,i)),'filled','k'); 
   line([0.1*i-0.02 0.1*i+0.02],[mean(GpkMax(:,i)) mean(GpkMax(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlabel('vt (cm/sec)');

subplot(2,2,2);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+0.1*i,22.5*squeeze(GpkWid(:,i)),'filled','k');
   line([0.1*i-0.02 0.1*i+0.02],[22.5*mean(GpkWid(:,i)) 22.5*mean(GpkWid(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
ylim([0 250]);
xlabel('vt (cm/sec)');

subplot(2,2,3);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+0.1*i,squeeze(RpkMax(:,i)),'filled','k'); 
   line([0.1*i-0.02 0.1*i+0.02],[mean(RpkMax(:,i)) mean(RpkMax(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak max');
xlabel('vt (cm/sec)');

subplot(2,2,4);
hold on;
for i=1:length(velThreshLC)
   scatter(zeros(numFlies,1)+0.1*i,22.5*squeeze(RpkWid(:,i)),'filled','k');
   line([0.1*i-0.02 0.1*i+0.02],[22.5*mean(RpkWid(:,i)) 22.5*mean(RpkWid(:,i))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
ylim([0 250]);
xlabel('vt (cm/sec)');


tStop = 5;
vRThresh = 0.075;
numPts = 17;

for flyID = 1:length(allFlyData)
    GPlotsBegR = zeros(100,18,9);
    GPlotsEndR = zeros(100,18,9);
    GPlotsBegL = zeros(100,18,9);
    GPlotsEndL = zeros(100,18,9);
    RPlotsBegR = zeros(100,18,9);
    RPlotsEndR = zeros(100,18,9);
    RPlotsBegL = zeros(100,18,9);
    RPlotsEndL = zeros(100,18,9);
    RitBegR = 1;
    RitBegL = 1;
    RitEndR = 1;
    RitEndL = 1;

    % Plot the rotational velocity vs the bump amplitude
    for darkID = 1:length(allFlyData{flyID}.dark)

        vRot = allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot;
        vF = diff(allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetForMatch);
        vZero = find(abs(vRot)<vRThresh);
        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);

        vZShort = [];
        for i=1:length(vZBeg)
            if vZEnd(i) - vZBeg(i) < tStop
                vZShort = horzcat(vZShort,[vZBeg(i):vZEnd(i)]);
            end
        end   
        vZero(vZShort) = [];

        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
%         vZEnd = vertcat(vZEnd,length(vZero));
%         vZBeg = vertcat(1,vZBeg);

    %     subplot(length(allFlyData{flyID}.dark),1,darkID)
    %     hold on;
    %     plot(vRot,'k');
    %     plot(vF*100,'m');
    %     scatter(vZero,vRot(vZero),'b');  
    %     scatter(vZero(vZEnd),vRot(vZero(vZEnd)),'r','filled');
    %     scatter(vZero(vZBeg),vRot(vZero(vZBeg)),'g','filled');

        if ~isempty(vZero)
            RSig = allFlyData{flyID}.dark{darkID}.RROIaveMax;
            GSig = allFlyData{flyID}.dark{darkID}.GROIaveMax;
            for i = 1:length(vZBeg)
                RSigNow = RSig(:,vZero(vZBeg(i)));
                GSigNow = GSig(:,vZero(vZBeg(i)));
%                 pkMax = max(RSigNow(1:8));
                pkMax = max(GSigNow(2:9));
%                 pkLoc = find(pkMax == RSigNow);
                pkLoc = find(pkMax == GSigNow);
                for j = 1:numPts
                    if vZero(vZBeg(i))-ceil(numPts/2)+j > 1 & vZero(vZBeg(i))-1 > 0 & vZero(vZBeg(i))-ceil(numPts/2)+j < length(RSig)
                        RSigNow(1:8) = circshift(RSig(1:8,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+3);
                        RSigNow(11:18) = circshift(RSig(11:18,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+3);
                        GSigNow(2:9) = circshift(GSig(2:9,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+3);
                        GSigNow(10:17) = circshift(GSig(10:17,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+3);

                        if vRot(vZero(vZBeg(i))-1) < 0
                            RPlotsBegL(RitBegL,:,j) = RSigNow;
                            GPlotsBegL(RitBegL,:,j) = GSigNow;
                        else
                            RPlotsBegR(RitBegR,:,j) = RSigNow;
                            GPlotsBegR(RitBegR,:,j) = GSigNow;
                        end
                    end
                end
                if vRot(vZero(vZBeg(i))-1) < 0
                    RitBegL = RitBegL + 1;
                else
                    RitBegR = RitBegR + 1;
                end
            end
            for i = 1:length(vZEnd)
                RSigNow = RSig(:,vZero(vZEnd(i)));
                GSigNow = GSig(:,vZero(vZEnd(i)));
%                 pkMax = max(RSigNow(1:8));
                pkMax = max(GSigNow(2:9));
%                 pkLoc = find(pkMax == RSigNow);
                pkLoc = find(pkMax == GSigNow);
                for j = 1:numPts
                    if vZero(vZEnd(i))-ceil(numPts/2)+j <= length(RSig) & vZero(vZEnd(i))-ceil(numPts/2)+j > 1 & vZero(vZEnd(i))+1 < length(vRot)
                        RSigNow(1:8) = circshift(RSig(1:8,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+3);
                        RSigNow(11:18) = circshift(RSig(11:18,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+3);
                        GSigNow(2:9) = circshift(GSig(2:9,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+3);
                        GSigNow(10:17) = circshift(GSig(10:17,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+3);

                        if vRot(vZero(vZEnd(i))+1) < 0
                            RPlotsEndL(RitEndL,:,j) = RSigNow;
                            GPlotsEndL(RitEndL,:,j) = GSigNow;
                        else
                            RPlotsEndR(RitEndR,:,j) = RSigNow;
                            GPlotsEndR(RitEndR,:,j) = GSigNow;
                        end
                    end

                end
                if vRot(vZero(vZEnd(i))+1) < 0 
                    RitEndL = RitEndL + 1;
                else
                    RitEndR = RitEndR + 1;
                end
            end
        end
    end

    RPlotsBegR(RitBegR:100,:,:) = [];
    RPlotsBegL(RitBegL:100,:,:) = [];
    RPlotsEndR(RitEndR:100,:,:) = [];
    RPlotsEndL(RitEndL:100,:,:) = [];
    GPlotsBegR(RitBegR:100,:,:) = [];
    GPlotsBegL(RitBegL:100,:,:) = [];
    GPlotsEndR(RitEndR:100,:,:) = [];
    GPlotsEndL(RitEndL:100,:,:) = [];

%     figure;
%     for i = 1:numPts
%         subplot(numPts,2,2*i-1);
%         hold on;
%         plot(quantile(squeeze(RPlotsBeg(:,:,i)),0.25,1),'r');
%         plot(quantile(squeeze(RPlotsBeg(:,:,i)),0.75,1),'r');
%         plot(mean(squeeze(RPlotsBeg(:,:,i)),1),'r','LineWidth',2);
%         plot(quantile(squeeze(GPlotsBeg(:,:,i)),0.25,1),'g');
%         plot(quantile(squeeze(GPlotsBeg(:,:,i)),0.75,1),'g');
%         plot(mean(squeeze(GPlotsBeg(:,:,i)),1),'g','LineWidth',2);
%         ylim([1 2]);
%         
%         subplot(numPts,2,2*i);
%         hold on;
%         plot(quantile(squeeze(RPlotsEnd(:,:,i)),0.25,1),'r');
%         plot(quantile(squeeze(RPlotsEnd(:,:,i)),0.75,1),'r');
%         plot(mean(squeeze(RPlotsEnd(:,:,i)),1),'r','LineWidth',2);
%         plot(quantile(squeeze(GPlotsEnd(:,:,i)),0.25,1),'g');
%         plot(quantile(squeeze(GPlotsEnd(:,:,i)),0.75,1),'g');
%         plot(mean(squeeze(GPlotsEnd(:,:,i)),1),'g','LineWidth',2);
%         ylim([1 2]);
%     end

    StartStopFig = figure;
    tStep = mean(diff(allFlyData{flyID}.dark{1}.positionDatMatch.OffsetRotMatch(:,1)));
    RMaxR = max(max(max(mean(squeeze(RPlotsBegR(:,[1:8 11:18],:)),1))),max(max(mean(squeeze(RPlotsEndR(:,[1:8 11:18],:)),1))));
    RMinR = min(min(min(mean(squeeze(RPlotsBegR(:,[1:8 11:18],:)),1))),min(min(mean(squeeze(RPlotsEndR(:,[1:8 11:18],:)),1))));
    RMaxL = max(max(max(mean(squeeze(RPlotsBegL(:,[1:8 11:18],:)),1))),max(max(mean(squeeze(RPlotsEndL(:,[1:8 11:18],:)),1))));
    RMinL = min(min(min(mean(squeeze(RPlotsBegL(:,[1:8 11:18],:)),1))),min(min(mean(squeeze(RPlotsEndL(:,[1:8 11:18],:)),1))));
    GMaxR = max(max(max(mean(squeeze(GPlotsBegR(:,2:17,:)),1))),max(max(mean(squeeze(GPlotsEndR(:,2:17,:)),1))));
    GMinR = min(min(min(mean(squeeze(GPlotsBegR(:,2:17,:)),1))),min(min(mean(squeeze(GPlotsEndR(:,2:17,:)),1))));
    GMaxL = max(max(max(mean(squeeze(GPlotsBegL(:,2:17,:)),1))),max(max(mean(squeeze(GPlotsEndL(:,2:17,:)),1))));
    GMinL = min(min(min(mean(squeeze(GPlotsBegL(:,2:17,:)),1))),min(min(mean(squeeze(GPlotsEndL(:,2:17,:)),1))));
    for i = 1:numPts
        RIm = mean(squeeze(RPlotsBegR(:,:,i)),1);
        RIm = (RIm-RMinR)./(RMaxR-RMinR);
        RIm(9:10) = 0;
        GIm = mean(squeeze(GPlotsBegR(:,:,i)),1);
        GIm = (GIm-GMinR)./(GMaxR-GMinR);
        GIm([1 18]) = 0;
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,4,4*i-3);
        imagesc(AllIm);
        caxis([0 1]);
        axis off;
        title(strcat('t_{Stop}=',num2str(tStep*(i-ceil(numPts/2))),' R turns'));
        
        RIm = mean(squeeze(RPlotsBegL(:,:,i)),1);
        RIm = (RIm-RMinL)./(RMaxL-RMinL);
        RIm(9:10) = 0;
        GIm = mean(squeeze(GPlotsBegL(:,:,i)),1);
        GIm = (GIm-GMinL)./(GMaxR-GMinL);
        GIm([1 18]) = 0;
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,4,4*i-2);
        imagesc(AllIm);
        caxis([0 1]);
        axis off;
        title(strcat('t_{Stop}=',num2str(tStep*(i-ceil(numPts/2))),' L turns'));
        
        RIm = mean(squeeze(RPlotsEndR(:,:,i)),1);
        RIm = (RIm-RMinR)./(RMaxR-RMinR);
        RIm(9:10) = 0;
        GIm = mean(squeeze(GPlotsEndR(:,:,i)),1);
        GIm = (GIm-GMinR)./(GMaxR-GMinR);
        GIm([1 18]) = 0;
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,4,4*i-1);
        imagesc(AllIm);
        axis off;
        title(strcat('t_{Start}=',num2str(tStep*(i-ceil(numPts/2))), 'R turns'));
        
        RIm = mean(squeeze(RPlotsEndL(:,:,i)),1);
        RIm = (RIm-RMinL)./(RMaxL-RMinL);
        RIm(9:10) = 0;
        GIm = mean(squeeze(GPlotsEndL(:,:,i)),1);
        GIm = (GIm-GMinL)./(GMaxL-GMinL);
        GIm([1 18]) = 0;
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,4,4*i);
        imagesc(AllIm);
        axis off;
        title(strcat('t_{Start}=',num2str(tStep*(i-ceil(numPts/2))), 'L turns'));
        
        colormap(brewermap(64, 'Blues'));
    end

    set(StartStopFig,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
    print(StartStopFig,strcat(allFlyData{flyID}.ID,'_StartStop_Galign'),'-dpdf');

end

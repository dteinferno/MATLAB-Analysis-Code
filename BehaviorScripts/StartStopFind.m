tStop = 5;
vRThresh = 0.075;
numPts = 17;

for flyID = 1:length(allFlyData)
    RPlotsBeg = zeros(100,16,9);
    RPlotsEnd = zeros(100,16,9);
    RitBeg = 1;
    RitEnd = 1;
    GPlotsBeg = zeros(100,16,9);
    GPlotsEnd = zeros(100,16,9);

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
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);

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
                pkMax = max(RSigNow);
                pkLoc = find(pkMax == RSigNow);
                for j = 1:numPts
                    if vZero(vZBeg(i))-ceil(numPts/2)+j > 1 & vZero(vZBeg(i))-1 > 0 & vZero(vZBeg(i))-ceil(numPts/2)+j < length(RSig)
                        RSigNow = circshift(RSig(:,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+8);
                        GSigNow = circshift(GSig(:,vZero(vZBeg(i))-ceil(numPts/2)+j),-pkLoc+8);

                        if vRot(vZero(vZBeg(i))-1) < 0
                            RPlotsBeg(RitBeg,:,j) = RSigNow;
                            GPlotsBeg(RitBeg,:,j) = GSigNow;
                        else
                            RPlotsBeg(RitBeg,:,j) = flipud(RSigNow);
                            GPlotsBeg(RitBeg,:,j) = flipud(GSigNow);
                        end
                    end
                end
                RitBeg = RitBeg + 1;
            end
            for i = 1:length(vZEnd)
                RSigNow = RSig(:,vZero(vZEnd(i)));
                GSigNow = GSig(:,vZero(vZEnd(i)));
                pkMax = max(RSigNow);
                pkLoc = find(pkMax == RSigNow);
                for j = 1:numPts
                    if vZero(vZEnd(i))-ceil(numPts/2)+j <= length(RSig) & vZero(vZEnd(i))-ceil(numPts/2)+j > 1 & vZero(vZEnd(i))+1 < length(vRot)
                        RSigNow = circshift(RSig(:,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+8);
                        GSigNow = circshift(GSig(:,vZero(vZEnd(i))-ceil(numPts/2)+j),-pkLoc+8);

                        if vRot(vZero(vZEnd(i))+1) < 0
                            RPlotsEnd(RitEnd,:,j) = RSigNow;
                            GPlotsEnd(RitEnd,:,j) = GSigNow;
                        else
                            RPlotsEnd(RitEnd,:,j) = flipud(RSigNow);
                            GPlotsEnd(RitEnd,:,j) = flipud(GSigNow);
                        end
                    end

                end
                RitEnd = RitEnd + 1;
            end
        end
    end

    RPlotsBeg(RitBeg:100,:,:) = [];
    GPlotsBeg(RitBeg:100,:,:) = [];
    RPlotsEnd(RitEnd:100,:,:) = [];
    GPlotsEnd(RitEnd:100,:,:) = [];

%     figure;
%     for i = 1:numPts
%         subplot(numPts,2,2*i-1);
%         hold on;
%         plot(quantile(squeeze(RPlotsBeg(:,:,i)),0.25,1),'r');
%         plot(quantile(squeeze(RPlotsBeg(:,:,i)),0.75,1),'r');
%         plot(median(squeeze(RPlotsBeg(:,:,i)),1),'r','LineWidth',2);
%         plot(quantile(squeeze(GPlotsBeg(:,:,i)),0.25,1),'g');
%         plot(quantile(squeeze(GPlotsBeg(:,:,i)),0.75,1),'g');
%         plot(median(squeeze(GPlotsBeg(:,:,i)),1),'g','LineWidth',2);
%         ylim([1 2]);
%         
%         subplot(numPts,2,2*i);
%         hold on;
%         plot(quantile(squeeze(RPlotsEnd(:,:,i)),0.25,1),'r');
%         plot(quantile(squeeze(RPlotsEnd(:,:,i)),0.75,1),'r');
%         plot(median(squeeze(RPlotsEnd(:,:,i)),1),'r','LineWidth',2);
%         plot(quantile(squeeze(GPlotsEnd(:,:,i)),0.25,1),'g');
%         plot(quantile(squeeze(GPlotsEnd(:,:,i)),0.75,1),'g');
%         plot(median(squeeze(GPlotsEnd(:,:,i)),1),'g','LineWidth',2);
%         ylim([1 2]);
%     end

    StartStopFig = figure;
    tStep = mean(diff(allFlyData{flyID}.dark{1}.positionDatMatch.OffsetRotMatch(:,1)));
    RMax = max(max(max(median(squeeze(RPlotsBeg(:,:,:)),1))),max(max(median(squeeze(RPlotsEnd(:,:,:)),1))));
    RMin = min(min(min(median(squeeze(RPlotsBeg(:,:,:)),1))),min(min(median(squeeze(RPlotsEnd(:,:,:)),1))));
    GMax = max(max(max(median(squeeze(GPlotsBeg(:,:,:)),1))),max(max(median(squeeze(GPlotsEnd(:,:,:)),1))));
    GMin = min(min(min(median(squeeze(GPlotsBeg(:,:,:)),1))),min(min(median(squeeze(GPlotsEnd(:,:,:)),1))));
    for i = 1:numPts
        RIm = median(squeeze(RPlotsBeg(:,:,i)),1);
        RIm = (RIm-RMin)./(RMax-RMin);
        GIm = median(squeeze(GPlotsBeg(:,:,i)),1);
        GIm = (GIm-GMin)./(GMax-GMin);
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,2,2*i-1);
        imagesc(AllIm);
        caxis([0 1]);
        axis off;
        title(strcat('t_{Stop}=',num2str(tStep*(i-ceil(numPts/2)))));
        
        RIm = median(squeeze(RPlotsEnd(:,:,i)),1);
        RIm = (RIm-RMin)./(RMax-RMin);
        GIm = median(squeeze(GPlotsEnd(:,:,i)),1);
        GIm = (GIm-GMin)./(GMax-GMin);
        AllIm = vertcat(RIm,GIm);
        subplot(numPts,2,2*i);
        imagesc(AllIm);
        axis off;
        title(strcat('t_{Start}=',num2str(tStep*(i-ceil(numPts/2)))));
        
        colormap(brewermap(64, 'Blues'));
    end

    set(StartStopFig,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
    print(StartStopFig,strcat(allFlyData{flyID}.ID,'_StartStop'),'-dpdf');

end

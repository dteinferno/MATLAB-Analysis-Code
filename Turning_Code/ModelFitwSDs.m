%% Get the difference between PVA direction for each fly
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=1:numFlies
    for velCut = 1:length(velThreshLC)       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
           circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),-1),...
           circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec),-1));
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
           circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),-1),...
           circshift(flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec),-1));
       
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
       angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
       angsraw = angsraw';
       PVAdiff(flyID,velCut) = -circ_mean(angsraw,mean(RturnsR,2))+circ_mean(angsraw,mean(RturnsG,2));

       PVAdiff(flyID,velCut) = min(PVAdiff(flyID,velCut),2*pi-PVAdiff(flyID,velCut));
       
    end

end

PVAdiff = PVAdiff;

PVAdiffplot = figure;
set(gcf,'Position',[100 100 560 420]);
subplot(5,1,[1:4]);
hold on;
% plot(180/pi*mean(PVAdiff,1),'Color','k','LineWidth',2);
plot(180/pi*(mean(PVAdiff,1)-std(PVAdiff,[],1)),'Color','k','LineWidth',2);
plot(180/pi*(mean(PVAdiff,1)+std(PVAdiff,[],1)),'Color','k','LineWidth',2);
set(gca,'FontSize',14);
ylabel('PVA difference (deg)');
xlim([1 6]);
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('v_r bins (deg/sec)');

ylim([-10 50]);
% ylim([-5 30]);

%%
PVAdiff = zeros(numFlies,length(velThreshLC));
for flyID=1:length(allFlyData)
    GallDat = [];
    RallDat = [];
    for trial = 1:length(allFlyData{flyID}.dark)
        GallDat = horzcat(GallDat,...
            reshape(allFlyData{flyID}.dark{trial}.GROIaveMax,...
            [1, size(allFlyData{flyID}.dark{trial}.GROIaveMax,1)*...
            size(allFlyData{flyID}.dark{trial}.GROIaveMax,2)]));
        RallDat = horzcat(RallDat,...
            reshape(allFlyData{flyID}.dark{trial}.RROIaveMax,...
            [1, size(allFlyData{flyID}.dark{trial}.RROIaveMax,1)*...
            size(allFlyData{flyID}.dark{trial}.RROIaveMax,2)]));
    end
    GSDNow = std(GallDat);
    RSDNow = std(RallDat);
    for velCut = 1:length(velThreshLC)     
       GturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkGalignDec));
       GturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkGalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkGalignDec));
       
       RturnsG = horzcat(XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsGAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsGAllDarkRalignDec));
       RturnsR = horzcat(XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignInc,...
           XCVals{flyID}.vel{velCut}.RturnsRAllDarkRalignDec,...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignInc),...
           flipud(XCVals{flyID}.vel{velCut}.LturnsRAllDarkRalignDec));

       if ~isempty(RturnsG) 
           % Calculate the PVA differences
           GpkMax(flyID, velCut) = mean(max(RturnsG,[],1)-min(RturnsG,[],1));
           GpkWids = zeros(size(RturnsG,2),1);
           
           for git = 1:size(RturnsG,2)
               XCNow = RturnsG(:,git);
               if max(XCNow) > mean(XCNow) + 0*GSDNow
                   overThresh = find((XCNow - min(XCNow))> 0.5*(max(XCNow)-min(XCNow)));
                   GpkWids(git) = length(overThresh)-1;
               end
           end
           GpkWid(flyID, velCut) = mean(GpkWids(find(GpkWids>0)));
               
           RpkMax(flyID, velCut) = mean(max(GturnsR,[],1)-min(GturnsR,[],1));
           RpkWids = zeros(size(GturnsR,2),1);
           for git = 1:size(GturnsR,2)
               XCNow = GturnsR(:,git);
               if max(XCNow) > mean(XCNow) + 0*RSDNow
                   overThresh = find((XCNow - min(XCNow))> 0.5*(max(XCNow)-min(XCNow)));
                   RpkWids(git) = length(overThresh)-1;
               end
           end
           RpkWid(flyID, velCut) = mean(RpkWids(find(RpkWids>0)));
       end
    end

end


GNorm = GpkMax(:, 1);
RNorm = RpkMax(:, 1);
for flyID = 1:5
    for velCut = 1:length(velThreshLC)
        GpkMax(flyID,velCut) = GpkMax(flyID,velCut)/mean(GNorm(1:5));
        RpkMax(flyID,velCut) = RpkMax(flyID,velCut)/mean(RNorm(1:5));
    end
end
for flyID = 6:10
    for velCut = 1:length(velThreshLC)
        GpkMax(flyID,velCut) = GpkMax(flyID,velCut)/mean(GNorm(6:10));
        RpkMax(flyID,velCut) = RpkMax(flyID,velCut)/mean(RNorm(6:10));
    end
end
        

PVAdiff = PVAdiff;

figure;
set(gcf,'Position',[100 100 560*2.5 420*2]);
subplot(2,2,1);
hold on;
GpkMaxMean = [];
GpkMaxSD = [];
for i=1:length(velThreshLC)
    PltNow = GpkMax(:,i);
    GpkMaxMean(i) = mean(PltNow(find(PltNow~=0)));
    GpkMaxSD(i) = std(PltNow(find(PltNow~=0)));
end
plot([15:30:30*length(velThreshLC)-15],GpkMaxMean+GpkMaxSD,'k');
plot([15:30:30*length(velThreshLC)-15],GpkMaxMean-GpkMaxSD,'k');
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 180]);
ylim([0.5 2.5]);
xlabel('vr bins (deg/sec)');

subplot(2,2,2);
hold on;
GpkWidMean = [];
GpkWidSD = [];
for i=1:length(velThreshLC)
   PltNow = GpkWid(:,i);
   GpkWidMean(i) = mean(PltNow(find(PltNow~=0)));
   GpkWidSD(i) = std(PltNow(find(PltNow~=0)));
end
plot([15:30:30*length(velThreshLC)-15],22.5*(GpkWidMean+GpkWidSD),'k');
plot([15:30:30*length(velThreshLC)-15],22.5*(GpkWidMean-GpkWidSD),'k');
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 180]);
xlabel('vr bins (deg/sec)');
ylim([0 180]);

subplot(2,2,3);
hold on;
RpkMaxMean = [];
RpkMaxSD = [];
for i=1:length(velThreshLC)
   PltNow = RpkMax(:,i);
   RpkMaxMean(i) = mean(PltNow(find(PltNow~=0)));
   RpkMaxSD(i) = std(PltNow(find(PltNow~=0)));
end
plot([15:30:30*length(velThreshLC)-15],RpkMaxMean+RpkMaxSD,'k');
plot([15:30:30*length(velThreshLC)-15],RpkMaxMean-RpkMaxSD,'k');
set(gca,'FontSize',14);
ylabel('Peak max');
xlim([0 180]);
ylim([0.5 2.5]);
xlabel('vr bins (deg/sec)');

subplot(2,2,4);
hold on;
RpkWidMean = [];
RpkWidSD = [];
for i=1:length(velThreshLC)
   PltNow = RpkWid(:,i);
   RpkWidMean(i) = mean(PltNow(find(PltNow~=0)));
   RpkWidSD(i) = std(PltNow(find(PltNow~=0)));
end
plot([15:30:30*length(velThreshLC)-15],22.5*(RpkWidMean+RpkWidSD),'k');
plot([15:30:30*length(velThreshLC)-15],22.5*(RpkWidMean-RpkWidSD),'k');
set(gca,'FontSize',14);
ylabel('Peak width (deg)');
xlim([0 180]);
xlabel('vr bins (deg/sec)');
ylim([0 180]);

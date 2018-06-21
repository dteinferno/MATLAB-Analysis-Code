%% Get the difference between peaks for each fly
PVAdiff = zeros(numFlies,length(velThreshLC));
PVAdiffMean = zeros(numFlies,length(velThreshLC));
PVAdiffSD = zeros(numFlies,length(velThreshLC));
Maxdiff = zeros(numFlies,length(velThreshLC));
MaxdiffMean = zeros(numFlies,length(velThreshLC));
MAxdiffSD = zeros(numFlies,length(velThreshLC));

num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

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

       if ~isempty(RturnsR)
           % Calculate the PVA differences for the mean signals
           PVAdiff(flyID,velCut) = -circ_mean(angsraw,mean(RturnsR,2))+circ_mean(angsraw,mean(RturnsG,2));
           PVAdiff(flyID,velCut) = min(PVAdiff(flyID,velCut),2*pi-PVAdiff(flyID,velCut));

           % Calculate the PVA differences for each signal and find the mean
           % and SD of these values
           PVAdiffTemp = zeros(1,size(RturnsR,2));
           for i = 1:size(RturnsR,2)
               PVAdiffTemp(i) = -circ_mean(angsraw,RturnsR(:,i))+circ_mean(angsraw,RturnsG(:,i));
               PVAdiffTemp(i) = min(PVAdiffTemp(i),2*pi-PVAdiffTemp(i));
           end
           PVAdiffMean(flyID,velCut) = mean(PVAdiffTemp);
           PVAdiffSD(flyID,velCut) = std(PVAdiffTemp);

           % Calculate the max differences for the mean signals
           RturnsRmean = mean(RturnsR,2);
           RturnsGmean = mean(RturnsG,2);
           Rmax = max(RturnsRmean);
           Gmax = max(RturnsGmean);
           Maxdiff(flyID,velCut) = -find(RturnsRmean==Rmax)+find(RturnsGmean==Gmax);
           Maxdiff(flyID,velCut) = min(Maxdiff(flyID,velCut),num_ROIs-Maxdiff(flyID,velCut));
           
           % Same as above for each signal
           MaxdiffTemp = zeros(1,size(RturnsR,2));
           for i = 1:size(RturnsR,2)
               Rmax = max(RturnsR(:,i));
               Gmax = max(RturnsG(:,i));
               MaxdiffTemp(i) = -find(RturnsR(:,i)==Rmax)+find(RturnsG(:,i)==Gmax);
               MaxdiffTemp(i) = min(MaxdiffTemp(i),num_ROIs-MaxdiffTemp(i));
           end
           MaxdiffMean(flyID,velCut) = mean(MaxdiffTemp);
           MaxdiffStd(flyID,velCut) = std(MaxdiffTemp);
       else
           flyID
           velCut
%            Maxdiff(flyID,velCut) = NaN;
%            MaxdiffMean(flyID,velCut) = NaN;
       end
       
       
    end

end

PVAdiff = PVAdiff;
PVAdiffMean = PVAdiffMean;
Maxdiff = Maxdiff;
MaxdiffMean = MaxdiffMean;

PVAdiffplot = figure;
set(gcf,'Position',[100 0 560*2 420*2]);
subplot(2,2,1);
hold on;
for plt = 1:6
    scatter(zeros(numFlies,1)+plt-1,180/pi*PVAdiff(:,plt),'filled','k');
    line([plt-1.25 plt-0.75],[180/pi*mean(PVAdiff(:,plt)) 180/pi*mean(PVAdiff(:,plt))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('PVA difference of mean (o)');
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('vRot bins (o/s)');
xlim([-0.5 5.5]);
ylim([-10 50]);

'PVA mean'
180/pi*mean(PVAdiff(:,plt))
180/pi*std(PVAdiff(:,plt))


subplot(2,2,2);
hold on;
for plt = 1:6
    scatter(zeros(numFlies,1)+plt-1,180/pi*PVAdiffMean(:,plt),'filled','k');
    line([plt-1.25 plt-0.75],[180/pi*mean(PVAdiffMean(:,plt)) 180/pi*mean(PVAdiffMean(:,plt))],'Color','k','LineWidth',2);
%     line([plt-1 plt-1],[180/pi*(PVAdiffMean(:,plt)-PVAdiffSD(:,plt)) 180/pi*(PVAdiffMean(:,plt)+PVAdiffSD(:,plt))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('mean PVA difference (o)');
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('vRot bins (o/sec)');
xlim([-0.5 5.5]);
ylim([-10 50]);

'PVA ind'
180/pi*mean(PVAdiffMean(:,plt))
180/pi*std(PVAdiffMean(:,plt))

subplot(2,2,3);
hold on;
for plt = 1:6
    scatter(zeros(numFlies,1)+plt-1,22.5*Maxdiff(:,plt),'filled','k');
    line([plt-1.25 plt-0.75],[22.5*mean(Maxdiff(:,plt)) 22.5*mean(Maxdiff(:,plt))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('peak difference of mean (o)');
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('vRot bins (o/s)');
xlim([-0.5 5.5]);
ylim([-10 50]);

'peak mean'
22.5*mean(Maxdiff(:,plt))
22.5*std(Maxdiff(:,plt))

subplot(2,2,4);
hold on;
for plt = 1:6
    scatter(zeros(numFlies,1)+plt-1,22.5*MaxdiffMean(:,plt),'filled','k');
    line([plt-1.25 plt-0.75],[22.5*mean(MaxdiffMean(:,plt)) 22.5*mean(MaxdiffMean(:,plt))],'Color','k','LineWidth',2);
%     line([plt-1 plt-1],[180/pi*(PVAdiffMean(:,plt)-PVAdiffSD(:,plt)) 180/pi*(PVAdiffMean(:,plt)+PVAdiffSD(:,plt))],'Color','k','LineWidth',2);
end
set(gca,'FontSize',14);
ylabel('mean peak difference (o)');
set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
xlabel('vRot bins (o/s)');
xlim([-0.5 5.5]);
ylim([-10 50]);

'peak ind'
22.5*mean(MaxdiffMean(:,plt))
22.5*std(MaxdiffMean(:,plt))

% set(PVAdiffplot,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
% print(PVAdiffplot,strcat('Results\','PVADiffnStrenOtherColor'),'-dpdf');

%% Plot examples for metrics
num_ROIs = 16;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';

RturnsRmean = mean(RturnsR,2);
RturnsRstd = std(RturnsR,[],2);
RturnsGmean = mean(RturnsG,2);
RturnsGstd = std(RturnsG,[],2);

exPlots = figure('units','normalized','outerposition',[0 0 1 0.25]);
subplot(1,4,1)
hold on;
plot(angsraw,(RturnsRmean-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsRmean+RturnsRstd-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsRmean-RturnsRstd-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsGmean-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'g');
plot(angsraw,(RturnsGmean+RturnsGstd-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'r');
plot(angsraw,(RturnsGmean-RturnsGstd-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'r');
line([circ_mean(angsraw,mean(RturnsR,2)) circ_mean(angsraw,mean(RturnsR,2))],...
    [0 1],'Color','r');
line([circ_mean(angsraw,mean(RturnsG,2)) circ_mean(angsraw,mean(RturnsG,2))],...
    [0 1],'Color','g');
xlim([-pi pi]);
ylim([-0.5 1.5]);

subplot(1,4,2)
hold on;
for i = 1:5%size(RturnsR,2);
    plot(angsraw,(RturnsR(:,i)-min(RturnsR(:,i)))./(max(RturnsR(:,i))-min(RturnsR(:,i))),'r');
    plot(angsraw,(RturnsG(:,i)-min(RturnsG(:,i)))./(max(RturnsG(:,i))-min(RturnsG(:,i))),'g');
    line([circ_mean(angsraw,mean(RturnsR(:,i),2)) circ_mean(angsraw,mean(RturnsR(:,i),2))],...
        [0 1],'Color','r');
    line([circ_mean(angsraw,mean(RturnsG(:,i),2)) circ_mean(angsraw,mean(RturnsG(:,i),2))],...
        [0 1],'Color','g');
end
xlim([-pi pi]);
ylim([-0.5 1.5]);

subplot(1,4,3)
hold on;
plot(angsraw,(RturnsRmean-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsRmean+RturnsRstd-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsRmean-RturnsRstd-min(RturnsRmean))./(max(RturnsRmean)-min(RturnsRmean)),'r');
plot(angsraw,(RturnsGmean-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'g');
plot(angsraw,(RturnsGmean+RturnsGstd-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'r');
plot(angsraw,(RturnsGmean-RturnsGstd-min(RturnsGmean))./(max(RturnsGmean)-min(RturnsGmean)),'r');
line(pi/8*[find(RturnsRmean==max(RturnsRmean)) find(RturnsRmean==max(RturnsRmean))]-pi,...
    [0 1],'Color','r');
line(pi/8*[find(RturnsGmean==max(RturnsGmean)) find(RturnsGmean==max(RturnsGmean))]-pi,...
    [0 1],'Color','g');
xlim([-pi pi]);
ylim([-0.5 1.5]);

subplot(1,4,4)
hold on;
for i = 1:5%size(RturnsR,2);
    plot(angsraw,(RturnsR(:,i)-min(RturnsR(:,i)))./(max(RturnsR(:,i))-min(RturnsR(:,i))),'r');
    plot(angsraw,(RturnsG(:,i)-min(RturnsG(:,i)))./(max(RturnsG(:,i))-min(RturnsG(:,i))),'g');
    line(pi/8*[find(RturnsR(:,i)==max(RturnsR(:,i))) find(RturnsR(:,i)==max(RturnsR(:,i)))]-pi,...
        [0 1],'Color','r');
    line(pi/8*[find(RturnsG(:,i)==max(RturnsG(:,i))) find(RturnsG(:,i)==max(RturnsG(:,i)))]-pi,...
        [0 1],'Color','g');
end
xlim([-pi pi]);
ylim([-0.5 1.5]);
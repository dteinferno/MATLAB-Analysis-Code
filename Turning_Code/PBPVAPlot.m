load('D:\Imaging\2Color\20160621\PVADiff_60D05-37F062.mat')
LWmT = LPVAdiff;
TLStrenWmT = LRPVAstren;
WLStrenWmT = LGPVAstren;
RWmT = RPVAdiff;
TRStrenWmT = RRPVAstren;
WRStrenWmT = RGPVAstren;
load('D:\Imaging\2Color\20160622\PVADiff_37F06-60D052.mat')
LTmW= LPVAdiff;
TLStrenTmW = LGPVAstren;
WLStrenTmW = LRPVAstren;
RTmW = RPVAdiff;
TRStrenTmW = RGPVAstren;
WRStrenTmW = RRPVAstren;


vRMin = 0;  
vRMax = 720;
vRSpan = 30;
plotRange = 150/vRSpan;

PVA = figure;
set(gcf,'Position',[100 100 560 420]);
subplot(5,1,[1:4]);
hold on;
for bnz = 1:(length(LWmT)-1)/2
    LNow = LWmT(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    LNow = -LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'r');
    RNow = RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'g');
    RNow = -RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
end
bnz = (length(LWmT)+1)/2;
LNow = LWmT(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
LNow = -LTmW(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'r');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz)))],'Color','r','LineWidth',3);
RNow = RWmT(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'g');
RNow = -RTmW(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz)))],'Color','r','LineWidth',3);

for bnz = (length(LWmT)+1)/2+1:length(LWmT)
    LNow = -LWmT(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'g');
    LNow = LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    RNow = RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'r');
    
%     scatter(zeros(6,1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,-LWmT(:,bnz),'r');
%     scatter(zeros(5,1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LTmW(:,bnz),'r');
%     scatter(zeros(6,1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,-RWmT(:,bnz),'g');
%     scatter(zeros(5,1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RTmW(:,bnz),'g');
end

% for flyID = 1:6
%     plot([(-plotRange*vRSpan):vRSpan:(plotRange*vRSpan)],-LWmT(flyID,:),'g');
%     plot([(-plotRange*vRSpan):vRSpan:(plotRange*vRSpan)],-RWmT(flyID,:),'r');
% end
% 
% for flyID = 1:5
%     plot([(-plotRange*vRSpan):vRSpan:(plotRange*vRSpan)],LTmW(flyID,:),'g','LineWidth',2);
%     plot([(-plotRange*vRSpan):vRSpan:(plotRange*vRSpan)],RTmW(flyID,:),'r','LineWidth',2);
% end
% 
for bnz = 1:(length(LWmT)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmT(:,bnz),-LTmW(:,bnz)),-RWmT(:,length(LWmT)+1-bnz)),RTmW(:,length(LWmT)+1-bnz));
%     allVal = vertcat(LWmT(:,bnz),-RWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    allVal = vertcat(vertcat(vertcat(RWmT(:,bnz),-RTmW(:,bnz))),-LWmT(:,length(LWmT)+1-bnz),LTmW(:,length(LWmT)+1-bnz));
%     allVal = vertcat(RWmT(:,bnz),-LWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);

subplot(5,1,5);
hold on;
% for bnz = 1:6
%     scatter(zeros(6,1)-(bnz-1)*vRSpan+plotRange*vRSpan,LStrenWmT(:,bnz),'r');
%     scatter(zeros(6,1)-(bnz-1)*vRSpan+plotRange*vRSpan,RStrenWmT(:,bnz),'r');
%     line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2 ...
%         -(bnz-1)*vRSpan+plotRange*vRSpan+vRSpan/2],...
%         [mean(vertcat(LStrenWmT(:,bnz),RStrenWmT(:,bnz))) ...
%         mean(vertcat(LStrenWmT(:,bnz),RStrenWmT(:,bnz)))],...
%         'Color','r');
%     scatter(zeros(5,1)-(bnz-1)*vRSpan+plotRange*vRSpan,LStrenTmW(:,bnz),'g');
%     scatter(zeros(5,1)-(bnz-1)*vRSpan+plotRange*vRSpan,RStrenTmW(:,bnz),'g');
%     line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2 ...
%         -(bnz-1)*vRSpan+plotRange*vRSpan+vRSpan/2],...
%         [mean(vertcat(LStrenTmW(:,bnz),RStrenTmW(:,bnz))) ...
%         mean(vertcat(LStrenTmW(:,bnz),RStrenTmW(:,bnz)))],...
%         'Color','g');
% end
% for bnz = 8:13
%     scatter(zeros(6,1)+(bnz-1)*vRSpan-plotRange*vRSpan,LStrenWmT(:,bnz),'r');
%     scatter(zeros(6,1)+(bnz-1)*vRSpan-plotRange*vRSpan,RStrenWmT(:,bnz),'r');
%     line([(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2 ...
%         (bnz-1)*vRSpan-plotRange*vRSpan+vRSpan/2],...
%         [mean(vertcat(LStrenWmT(:,bnz),RStrenWmT(:,bnz))) ...
%         mean(vertcat(LStrenWmT(:,bnz),RStrenWmT(:,bnz)))],...
%         'Color','r');
%     scatter(zeros(5,1)+(bnz-1)*vRSpan-plotRange*vRSpan,LStrenTmW(:,bnz),'g');
%     scatter(zeros(5,1)+(bnz-1)*vRSpan-plotRange*vRSpan,RStrenTmW(:,bnz),'g');
%     line([(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2 ...
%         (bnz-1)*vRSpan-plotRange*vRSpan+vRSpan/2],...
%         [mean(vertcat(LStrenTmW(:,bnz),RStrenTmW(:,bnz))) ...
%         mean(vertcat(LStrenTmW(:,bnz),RStrenTmW(:,bnz)))],...
%         'Color','g');
% end
% ylim([0 0.15]);
allValT = vertcat(vertcat(vertcat(vertcat(vertcat(vertcat(vertcat(...
    fliplr(TLStrenWmT(:,1:(length(LWmT)-1)/2)),...
    fliplr(TRStrenWmT(:,1:(length(LWmT)-1)/2))),...
    fliplr(TLStrenTmW(:,1:(length(LWmT)-1)/2))),...
    fliplr(TRStrenTmW(:,1:(length(LWmT)-1)/2))),...
    TLStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    TRStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    TLStrenTmW(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    TRStrenTmW(:,(length(LWmT)+1)/2+1:length(LWmT)));
allValTZero = vertcat(vertcat(vertcat(...
    TLStrenWmT(:,(length(LWmT)+1)/2),...
    TRStrenWmT(:,(length(LWmT)+1)/2)),...
    TLStrenTmW(:,(length(LWmT)+1)/2)),...
    TRStrenTmW(:,(length(LWmT)+1)/2));
% allValT(find(allValT==0)) = [];
allValW = vertcat(vertcat(vertcat(vertcat(vertcat(vertcat(vertcat(...
    fliplr(WLStrenWmT(:,1:(length(LWmT)-1)/2)),...
    fliplr(WRStrenWmT(:,1:(length(LWmT)-1)/2))),...
    fliplr(WLStrenTmW(:,1:(length(LWmT)-1)/2))),...
    fliplr(WRStrenTmW(:,1:(length(LWmT)-1)/2))),...
    WLStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    WRStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    WLStrenTmW(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    WRStrenTmW(:,(length(LWmT)+1)/2+1:length(LWmT)));
allValWZero = vertcat(vertcat(vertcat(...
    WLStrenWmT(:,(length(LWmT)+1)/2),...
    WRStrenWmT(:,(length(LWmT)+1)/2)),...
    WLStrenTmW(:,(length(LWmT)+1)/2)),...
    WRStrenTmW(:,(length(LWmT)+1)/2));
% allValW(find(allValW==0)) = [];
imagesc(vertcat(horzcat(mean(allValTZero),mean(allValT,1)),...
    horzcat(mean(allValWZero),mean(allValW,1))));
colormap(brewermap(64, 'Greys'));
xlim([0.5 (length(LWmT)-1)/2+1.5]);
caxis([0 0.05]);
colorbar;

% set(PVA,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
% print(PVA,strcat('../Results/PBPVASummary'),'-dpdf');

%% As above, but each point colored by PVA
% Load the data
load('D:\Imaging\2Color\20160621\PVADiff_60D05-37F062.mat')
LWmT = LPVAdiff;
TLStrenWmT = LRPVAstren;
WLStrenWmT = LGPVAstren;
RWmT = RPVAdiff;
TRStrenWmT = RRPVAstren;
WRStrenWmT = RGPVAstren;
load('D:\Imaging\2Color\20160622\PVADiff_37F06-60D052.mat')
LTmW= LPVAdiff;
TLStrenTmW = LGPVAstren;
WLStrenTmW = LRPVAstren;
RTmW = RPVAdiff;
TRStrenTmW = RGPVAstren;
WRStrenTmW = RRPVAstren;

% Specify the parameters used for sorting
vRMin = 0;  
vRMax = 720;
vRSpan = 30;
plotRange = 150/vRSpan;

% Specify the max PVA strength
PVAStren = 0.05;

% Make the figure
PVA = figure;
set(gcf,'Position',[100 100 560 420*2]);
subplot(2,1,1);
title('37F06');
hold on;
for bnz = 1:(length(LWmT)-1)/2
    LNow = LWmT(:,bnz);
    SNow = min(TLStrenWmT(:,bnz),WLStrenWmT(:,bnz));
    SNow(find(LNow==0)) = [];
    LNow(find(LNow==0)) = [];
    SCol = zeros(length(LNow),3);
    for it = 1:length(LNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/10,LNow,50,SCol,'filled');
    LNow = -LTmW(:,bnz);
    SNow = min(TLStrenTmW(:,bnz),WLStrenTmW(:,bnz));
    SNow(find(LNow==0)) = [];
    LNow(find(LNow==0)) = [];
    SCol = zeros(length(LNow),3);
    for it = 1:length(LNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/10,LNow,50,SCol,'filled');
    RNow = RWmT(:,bnz);
    SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
    SNow(find(RNow==0)) = [];
    RNow(find(RNow==0)) = [];
    SCol = zeros(length(RNow),3);
    for it = 1:length(RNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan+vRSpan/10,RNow,50,SCol,'LineWidth',1);
    RNow = -RTmW(:,bnz);
    SNow = min(TRStrenTmW(:,bnz),WRStrenTmW(:,bnz));
    SNow(find(RNow==0)) = [];
    RNow(find(RNow==0)) = [];
    SCol = zeros(length(RNow),3);
    for it = 1:length(RNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan+vRSpan/10,RNow,50,SCol,'LineWidth',1);
end
bnz = (length(LWmT)+1)/2;
LNow = LWmT(:,bnz);
SNow = min(TLStrenWmT(:,bnz),WLStrenWmT(:,bnz));
SNow(find(LNow==0)) = [];
LNow(find(LNow==0)) = [];
SCol = zeros(length(LNow),3);
for it = 1:length(LNow)
    SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
end
scatter(zeros(length(LNow),1)-vRSpan/10,LNow,50,SCol,'filled');
LNow = -LTmW(:,bnz);
SNow = min(TLStrenTmW(:,bnz),WLStrenTmW(:,bnz));
SNow(find(LNow==0)) = [];
LNow(find(LNow==0)) = [];
SCol = zeros(length(LNow),3);
for it = 1:length(LNow)
    SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
end
scatter(zeros(length(LNow),1)-vRSpan/10,LNow,50,SCol,'filled');
line([-vRSpan/4-vRSpan/10 vRSpan/4-vRSpan/10],[mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz)))],'Color','k','LineWidth',3);
RNow = RWmT(:,bnz);
SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
SNow(find(RNow==0)) = [];
RNow(find(RNow==0)) = [];
SCol = zeros(length(RNow),3);
for it = 1:length(RNow)
    SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
end
scatter(zeros(length(RNow),1)+vRSpan/10,RNow,50,SCol,'LineWidth',1);
RNow = -RTmW(:,bnz);
SNow = min(TRStrenTmW(:,bnz),WRStrenTmW(:,bnz));
SNow(find(RNow==0)) = [];
RNow(find(RNow==0)) = [];
SCol = zeros(length(RNow),3);
for it = 1:length(RNow)
    SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
end
scatter(zeros(length(RNow),1)+vRSpan/10,RNow,50,SCol,'LineWidth',1);
line([-vRSpan/4+vRSpan/10 vRSpan/4+vRSpan/10],[mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz)))],'Color','k','LineWidth',3);

for bnz = (length(LWmT)+1)/2+1:length(LWmT)
    LNow = -LWmT(:,bnz);
    SNow = min(TLStrenWmT(:,bnz),WLStrenWmT(:,bnz));
    SNow(find(LNow==0)) = [];
    LNow(find(LNow==0)) = [];
    SCol = zeros(length(LNow),3);
    for it = 1:length(LNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan+vRSpan/10,LNow,50,SCol,'LineWidth',1);
    LNow = LTmW(:,bnz);
    SNow = min(TLStrenTmW(:,bnz),WLStrenTmW(:,bnz));
    SNow(find(LNow==0)) = [];
    LNow(find(LNow==0)) = [];
    SCol = zeros(length(LNow),3);
    for it = 1:length(LNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan+vRSpan/10,LNow,50,SCol,'LineWidth',1);
    RNow = -RWmT(:,bnz);
    SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
    SNow(find(RNow==0)) = [];
    RNow(find(RNow==0)) = [];
    SCol = zeros(length(RNow),3);
    for it = 1:length(RNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/10,RNow,50,SCol,'filled');
    RNow = RTmW(:,bnz);
    SNow = min(TRStrenTmW(:,bnz),WRStrenTmW(:,bnz));
    SNow(find(RNow==0)) = [];
    RNow(find(RNow==0)) = [];
    SCol = zeros(length(RNow),3);
    for it = 1:length(RNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/10,RNow,50,SCol,'filled');
end

for bnz = 1:(length(LWmT)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmT(:,bnz),-LTmW(:,bnz)),-RWmT(:,length(LWmT)+1-bnz)),RTmW(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','k','LineWidth',3);
    allVal = vertcat(vertcat(vertcat(RWmT(:,bnz),-RTmW(:,bnz))),-LWmT(:,length(LWmT)+1-bnz),LTmW(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','k','LineWidth',3);
end
line([-vRSpan/4-vRSpan/10 vRSpan*plotRange+vRSpan/2],[0 0],'Color','k');
xlim([-vRSpan/4-vRSpan/10 165]);
ylim([-1.5 3]);

ylabel('PB offset (# of glom.)');

% set(PVA,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
% print(PVA,strcat('../Results/PBPVASummary'),'-dpdf');
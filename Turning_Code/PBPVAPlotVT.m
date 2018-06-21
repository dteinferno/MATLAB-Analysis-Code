load('D:\Imaging\2Color\20160903\PVADiff_60D05-VT81352.mat')
LWmT = LPVAdiff;
TLStrenWmT = LRPVAstren;
WLStrenWmT = LGPVAstren;
RWmT = RPVAdiff;
TRStrenWmT = RRPVAstren;
WRStrenWmT = RGPVAstren;

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
    RNow = RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
end
bnz = (length(LWmT)+1)/2;
LNow = LWmT(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(LWmT(:,bnz)) mean(LWmT(:,bnz))],'Color','g','LineWidth',3);
RNow = RWmT(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(RWmT(:,bnz)) mean(RWmT(:,bnz))],'Color','r','LineWidth',3);

for bnz = (length(LWmT)+1)/2+1:length(LWmT)
    LNow = -LWmT(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
end

for bnz = 1:(length(LWmT)-1)/2
    allVal = vertcat(LWmT(:,bnz),-RWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    allVal = vertcat(RWmT(:,bnz),-LWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);
subplot(5,1,5);
hold on;

allValT = vertcat(vertcat(vertcat(...
    fliplr(TLStrenWmT(:,1:(length(LWmT)-1)/2)),...
    fliplr(TRStrenWmT(:,1:(length(LWmT)-1)/2))),...
    TLStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    TRStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT)));
allValTZero = vertcat(...
    TLStrenWmT(:,(length(LWmT)+1)/2),...
    TRStrenWmT(:,(length(LWmT)+1)/2));
% allValT(find(allValT==0)) = [];
allValW = vertcat(vertcat(vertcat(...
    fliplr(WLStrenWmT(:,1:(length(LWmT)-1)/2)),...
    fliplr(WRStrenWmT(:,1:(length(LWmT)-1)/2))),...
    WLStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT))),...
    WRStrenWmT(:,(length(LWmT)+1)/2+1:length(LWmT)));
allValWZero = vertcat(...
    WLStrenWmT(:,(length(LWmT)+1)/2),...
    WRStrenWmT(:,(length(LWmT)+1)/2));
% allValW(find(allValW==0)) = [];
imagesc(vertcat(horzcat(mean(allValTZero),mean(allValT,1)),...
    horzcat(mean(allValWZero),mean(allValW,1))));
colormap(brewermap(64, 'Greys'));
xlim([0.5 (length(LWmT)-1)/2+1.5]);
caxis([0 0.05]);
colorbar;

set(PVA,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(PVA,strcat('../Results/PBPVASummary_VT8135'),'-dpdf');

%% As above, but each point colored by PVA
load('D:\Imaging\2Color\20160903\PVADiff_60D05-VT81352.mat')
LWmT = LPVAdiff;
TLStrenWmT = LRPVAstren;
WLStrenWmT = LGPVAstren;
RWmT = RPVAdiff;
TRStrenWmT = RRPVAstren;
WRStrenWmT = RGPVAstren;

vRMin = 0;  
vRMax = 720;
vRSpan = 30;
plotRange = 150/vRSpan;

% PVA = figure;
% set(gcf,'Position',[100 100 560 420]);
subplot(2,1,2);
hold on;
title('8135');
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
    RNow = RWmT(:,bnz);
    SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
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
% line([-vRSpan/4-vRSpan/10 vRSpan/4-vRSpan/10],[mean(LWmT(:,bnz)) mean(LWmT(:,bnz))],'Color','k','LineWidth',3);
RNow = RWmT(:,bnz);
SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
SNow(find(RNow==0)) = [];
RNow(find(RNow==0)) = [];
SCol = zeros(length(RNow),3);
for it = 1:length(RNow)
    SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
end
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)+vRSpan/10,RNow,50,SCol,'LineWidth',1);
% line([-vRSpan/4+vRSpan/10 vRSpan/4+vRSpan/10],[mean(RWmT(:,bnz)) mean(RWmT(:,bnz))],'Color','k','LineWidth',3);

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
    RNow = -RWmT(:,bnz);
    SNow = min(TRStrenWmT(:,bnz),WRStrenWmT(:,bnz));
    SNow(find(RNow==0)) = [];
    RNow(find(RNow==0)) = [];
    SCol = zeros(length(RNow),3);
    for it = 1:length(RNow)
        SCol(it,:) = (PVAStren-SNow(it))/PVAStren;
    end
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/10,RNow,50,SCol,'filled');
end

% for bnz = 1:(length(LWmT)-1)/2
%     allVal = vertcat(LWmT(:,bnz),-RWmT(:,length(LWmT)+1-bnz));
%     allVal(find(allVal==0)) = [];
%     line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/10],...
%         [mean(allVal) mean(allVal)], 'Color','k','LineWidth',3);
%     allVal = vertcat(RWmT(:,bnz),-LWmT(:,length(LWmT)+1-bnz));
%     allVal(find(allVal==0)) = [];
%     line([-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan+vRSpan/10],...
%         [mean(allVal) mean(allVal)], 'Color','k','LineWidth',3);
% end
line([-vRSpan/4-vRSpan/10 vRSpan*plotRange+vRSpan/2],[0 0],'Color','k');
xlim([-vRSpan/4-vRSpan/10 165]);
ylim([-1.5 3]);

ylabel('PB offset (# of glom.)');
xlabel('vRot (o/s)');

set(PVA,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(PVA,strcat('../Results/PBPVASummary_StrengthCoded'),'-dpdf');
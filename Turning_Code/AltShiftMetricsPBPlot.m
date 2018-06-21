load('D:\Imaging\2Color\20160621\PVADiff_60D05-37F062.mat')
LWmT = LPVAdiff;
LWmTMean = LPVAdiffMean;
LWmTPeak = LPeakdiff;
LWmTPeakMean = LPeakdiffMean;
RWmT = RPVAdiff;
RWmTMean = RPVAdiffMean;
RWmTPeak = RPeakdiff;
RWmTPeakMean = RPeakdiffMean;

load('D:\Imaging\2Color\20160622\PVADiff_37F06-60D052.mat')
LTmW= LPVAdiff;
LTmWMean = LPVAdiffMean;
LTmWPeak = LPeakdiff;
LTmWPeakMean = LPeakdiffMean;
RTmW = RPVAdiff;
RTmWMean = RPVAdiffMean;
RTmWPeak = RPeakdiff;
RTmWPeakMean = RPeakdiffMean;



vRMin = 0;  
vRMax = 720;
vRSpan = 30;
plotRange = 150/vRSpan;

PVA = figure;
set(gcf,'Position',[100 100 560 420]);

subplot(2,2,1);
hold on;
for bnz = 1:(length(LWmT)-1)/2
    LNow = LWmT(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    LNow = -LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    RNow = RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
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
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmT(:,bnz),-LTmW(:,bnz)))],'Color','g','LineWidth',3);
RNow = RWmT(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
RNow = -RTmW(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmT(:,bnz),-RTmW(:,bnz)))],'Color','r','LineWidth',3);

for bnz = (length(LWmT)+1)/2+1:length(LWmT)
    LNow = -LWmT(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    LNow = LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmT(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    RNow = RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    

end

for bnz = 1:(length(LWmT)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmT(:,bnz),-LTmW(:,bnz)),-RWmT(:,length(LWmT)+1-bnz)),RTmW(:,length(LWmT)+1-bnz));
%     allVal = vertcat(LWmT(:,bnz),-RWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    if bnz == 2
        'PVA mean'
        mean(allVal)
        std(allVal)
    end
    allVal = vertcat(vertcat(vertcat(RWmT(:,bnz),-RTmW(:,bnz))),-LWmT(:,length(LWmT)+1-bnz),LTmW(:,length(LWmT)+1-bnz));
%     allVal = vertcat(RWmT(:,bnz),-LWmT(:,length(LWmT)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
    if bnz == 2
        mean(allVal)
        std(allVal)
    end
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);
title('PVA diff of mean');

subplot(2,2,2);
hold on;
for bnz = 1:(length(LWmTMean)-1)/2
    LNow = LWmTMean(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    LNow = -LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    RNow = RWmTMean(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
    RNow = -RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
end
bnz = (length(LWmTMean)+1)/2;
LNow = LWmTMean(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
LNow = -LTmW(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(vertcat(LWmTMean(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmTMean(:,bnz),-LTmW(:,bnz)))],'Color','g','LineWidth',3);
RNow = RWmTMean(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
RNow = -RTmW(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(vertcat(RWmTMean(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmTMean(:,bnz),-RTmW(:,bnz)))],'Color','r','LineWidth',3);

for bnz = (length(LWmTMean)+1)/2+1:length(LWmTMean)
    LNow = -LWmTMean(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    LNow = LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmTMean(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    RNow = RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    

end

for bnz = 1:(length(LWmTMean)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmTMean(:,bnz),-LTmW(:,bnz)),-RWmTMean(:,length(LWmTMean)+1-bnz)),RTmW(:,length(LWmTMean)+1-bnz));
%     allVal = vertcat(LWmTMean(:,bnz),-RWmTMean(:,length(LWmTMean)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    if bnz == 2
        'PVA ind'
        mean(allVal)
        std(allVal)
    end
    allVal = vertcat(vertcat(vertcat(RWmTMean(:,bnz),-RTmW(:,bnz))),-LWmTMean(:,length(LWmTMean)+1-bnz),LTmW(:,length(LWmTMean)+1-bnz));
%     allVal = vertcat(RWmTMean(:,bnz),-LWmTMean(:,length(LWmTMean)+1-bnz));
    allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
    if bnz == 2
        mean(allVal)
        std(allVal)
    end
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);
title('mean of PVA diffs');

subplot(2,2,3);
hold on;
for bnz = 1:(length(LWmTPeak)-1)/2
    LNow = LWmTPeak(:,bnz);
%     LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    LNow = -LTmWPeak(:,bnz);
%     LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    RNow = RWmTPeak(:,bnz);
%     RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
    RNow = -RTmWPeak(:,bnz);
%     RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
end
bnz = (length(LWmTPeak)+1)/2;
LNow = LWmTPeak(:,bnz);
% LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
LNow = -LTmWPeak(:,bnz);
% LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(vertcat(LWmTPeak(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmTPeak(:,bnz),-LTmW(:,bnz)))],'Color','g','LineWidth',3);
RNow = RWmTPeak(:,bnz);
% RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
RNow = -RTmWPeak(:,bnz);
% RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(vertcat(RWmTPeak(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmTPeak(:,bnz),-RTmW(:,bnz)))],'Color','r','LineWidth',3);

for bnz = (length(LWmTPeak)+1)/2+1:length(LWmTPeak)
    LNow = -LWmTPeak(:,bnz);
%     LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    LNow = LTmWPeak(:,bnz);
%     LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmTPeak(:,bnz);
%     RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    RNow = RTmWPeak(:,bnz);
%     RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
end

for bnz = 1:(length(LWmTPeak)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmTPeak(:,bnz),-LTmW(:,bnz)),-RWmTPeak(:,length(LWmTPeak)+1-bnz)),RTmWPeak(:,length(LWmTPeak)+1-bnz));
%     allVal = vertcat(LWmTPeak(:,bnz),-RWmTPeak(:,length(LWmTPeak)+1-bnz));
%     allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    if bnz == 2
        'Peak mean'
        mean(allVal)
        std(allVal)
    end
    allVal = vertcat(vertcat(vertcat(RWmTPeak(:,bnz),-RTmW(:,bnz))),-LWmTPeak(:,length(LWmTPeak)+1-bnz),LTmWPeak(:,length(LWmTPeak)+1-bnz));
%     allVal = vertcat(RWmTPeak(:,bnz),-LWmTPeak(:,length(LWmTPeak)+1-bnz));
%     allVal(find(allVal==0)) = [];
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
    if bnz == 2
        mean(allVal)
        std(allVal)
    end
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);
title('peak diff of mean');

subplot(2,2,4);
hold on;
for bnz = 1:(length(LWmTPeakMean)-1)/2
    LNow = LWmTPeakMean(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    LNow = -LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'g');
    RNow = RWmTPeakMean(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
    RNow = -RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'r');
end
bnz = (length(LWmTPeakMean)+1)/2;
LNow = LWmTPeakMean(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
LNow = -LTmW(:,bnz);
LNow(find(LNow==0)) = [];
scatter(zeros(length(LNow),1)-vRSpan/2-vRSpan/10,LNow,'g');
line([-3*vRSpan/4-vRSpan/10 -vRSpan/4-vRSpan/10],[mean(vertcat(LWmTPeakMean(:,bnz),-LTmW(:,bnz))) mean(vertcat(LWmTPeakMean(:,bnz),-LTmW(:,bnz)))],'Color','g','LineWidth',3);
RNow = RWmTPeakMean(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
RNow = -RTmW(:,bnz);
RNow(find(RNow==0)) = [];
scatter(zeros(length(RNow),1)-vRSpan/2+vRSpan/10,RNow,'r');
line([-3*vRSpan/4+vRSpan/10 -vRSpan/4+vRSpan/10],[mean(vertcat(RWmTPeakMean(:,bnz),-RTmW(:,bnz))) mean(vertcat(RWmTPeakMean(:,bnz),-RTmW(:,bnz)))],'Color','r','LineWidth',3);

for bnz = (length(LWmTPeakMean)+1)/2+1:length(LWmTPeakMean)
    LNow = -LWmTPeakMean(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    LNow = LTmW(:,bnz);
    LNow(find(LNow==0)) = [];
    scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'r');
    RNow = -RWmTPeakMean(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    RNow = RTmW(:,bnz);
    RNow(find(RNow==0)) = [];
    scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'g');
    

end

for bnz = 1:(length(LWmTPeakMean)-1)/2
    allVal = vertcat(vertcat(vertcat(LWmTPeakMean(:,bnz),-LTmW(:,bnz)),-RWmTPeakMean(:,length(LWmTPeakMean)+1-bnz)),RTmW(:,length(LWmTPeakMean)+1-bnz));
%     allVal = vertcat(LWmTPeakMean(:,bnz),-RWmTPeakMean(:,length(LWmTPeakMean)+1-bnz));
    allVal(find(allVal==0)) = [];
    if bnz == 2
        'Peak ind'
        mean(allVal)
        std(allVal)
    end
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4-vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2-vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','g','LineWidth',3);
    allVal = vertcat(vertcat(vertcat(RWmTPeakMean(:,bnz),-RTmW(:,bnz))),-LWmTPeakMean(:,length(LWmTPeakMean)+1-bnz),LTmW(:,length(LWmTPeakMean)+1-bnz));
%     allVal = vertcat(RWmTPeakMean(:,bnz),-LWmTPeakMean(:,length(LWmTPeakMean)+1-bnz));
    allVal(find(allVal==0)) = [];
    if bnz == 2
        mean(allVal)
        std(allVal)
    end
    line([-(bnz-1)*vRSpan-vRSpan/2+plotRange*vRSpan-vRSpan/4+vRSpan/10,-(bnz-1)*vRSpan+vRSpan/4+plotRange*vRSpan-vRSpan/2+vRSpan/10],...
        [mean(allVal) mean(allVal)], 'Color','r','LineWidth',3);
end
line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
xlim([-3*vRSpan/4-vRSpan/10 150]);
ylim([-1.5 2.5]);
title('mean of peak diffs');


cd('D:\Imaging\2Color\Results');
load('NormDat');
figure;
hold on;
scatter(zeros(length(GNorm),1)+1,GNorm,'g','o','filled');
scatter(zeros(length(RNorm),1)+2,RNorm,'r','d','filled');
scatter(zeros(length(RNorm10(6:10)),1)+3,RNorm10(6:10),'r','o','filled');
scatter(zeros(length(GNorm10(6:10)),1)+4,GNorm10(6:10),'g','d','filled');
xlim([0 5]);

for i=1:5
    line([1 2],[GNorm(i) RNorm(i)]);
    line([3 4],[RNorm10(i+5) GNorm10(i+5)]);
end

set(gca,'XTick',[1:4]);
set(gca,'XTickLabel',{'wedges','tiles','wedges','tiles'});

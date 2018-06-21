RvG = figure;
subplot(3,2,1);
imagesc(squeeze(mean(RstackFull,3)))
caxis([0 1200]);
axis equal;
colorbar;
subplot(3,2,2);
imagesc(squeeze(mean(GstackFull,3)))
caxis([0 65]);
axis equal;
colorbar;
subplot(3,2,3);
imagesc(squeeze(mean(RstackFlat,3)))
caxis([0 300]);
axis equal;
colorbar;
subplot(3,2,4);
imagesc(squeeze(mean(GstackFlat,3)))
caxis([0 15]);
axis equal;
colorbar;
subplot(3,2,5);
imagesc(squeeze(mean(ROIs,1)))
axis equal;
colorbar;

blues = colormap(brewermap(64, 'Blues'));

% greens = blues;
% greens(:,1) = 1-blues(:,3);
% greens(:,2) = 1-blues(:,1);
% greens(:,3) = 1-blues(:,3);
% colormap(greens);

magentas(:,1) = 1-blues(:,1);
magentas(:,2) = 1-blues(:,3);
magentas(:,3) = 1-blues(:,1);
colormap(magentas);

set(RvG,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5  11]);
print(RvG,'RvsG_magentas','-dpdf');


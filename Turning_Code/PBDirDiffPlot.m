% Look at the red channel
LdirDiffR = zeros(length(allFlyData),5);
RdirDiffR = zeros(length(allFlyData),5);
LdirDiffG = zeros(length(allFlyData),5);
RdirDiffG = zeros(length(allFlyData),5);

for flyID = 1:length(allFlyData)
    
    for vRng = 1:5
        
        LdirDiffR(flyID,vRng) = mean(LPB.dark.dataR{flyID}.CW{vRng}) - mean(LPB.dark.dataR{flyID}.CCW{vRng});
        LdirDiffG(flyID,vRng) = mean(LPB.dark.dataG{flyID}.CW{vRng}) - mean(LPB.dark.dataG{flyID}.CCW{vRng});
        
        RdirDiffG(flyID,vRng) = mean(RPB.dark.dataG{flyID}.CW{vRng}) - mean(RPB.dark.dataG{flyID}.CCW{vRng});
        RdirDiffR(flyID,vRng) = mean(RPB.dark.dataR{flyID}.CW{vRng}) - mean(RPB.dark.dataR{flyID}.CCW{vRng});
        
    end
    
end


figure;
for vRng = 1:5
    subplot(1,2,1);
    hold on;
    scatter(zeros(length(allFlyData),1)+vRng,LdirDiffR(:,vRng),'r');
    scatter(zeros(length(allFlyData),1)+vRng,LdirDiffG(:,vRng),'g');
    ylim([-0.4 0.4]);
    line([1 5],[0 0],'Color','k');
    
    subplot(1,2,2);
    hold on;
    scatter(zeros(length(allFlyData),1)+vRng,RdirDiffR(:,vRng),'r');
    scatter(zeros(length(allFlyData),1)+vRng,RdirDiffG(:,vRng),'g');
    ylim([-0.4 0.4]);
    line([1 5],[0 0],'Color','k');
end

% save(strcat('DirDiff_',greenLine,'-',redLine,'.mat'),'LdirDiffR','RdirDiffR','LdirDiffG','RdirDiffG');

%%
load('D:\Imaging\2Color\20160621\DirDiff_60D05-37F06.mat');
WGL = LdirDiffG;
WGR = RdirDiffG;
TRL = LdirDiffR;
TRR = RdirDiffR;

load('D:\Imaging\2Color\20160622\DirDiff_37F06-60D05.mat');
WRL = LdirDiffR;
WRR = RdirDiffR;
TGL = LdirDiffG;
TGR = RdirDiffG;

figure;
subplot(2,2,1)
hold on;
for vRng = 1:5
    scatter(zeros(6,1)+vRng,TRL(:,vRng),'k','filled');
    scatter(zeros(5,1)+vRng,TGL(:,vRng),'k','filled');
    TLAll = vertcat(TRL(:,vRng),TGL(:,vRng));
    line([vRng-1/2 vRng+1/2],[mean(TLAll(~isnan(TLAll))) mean(TLAll(~isnan(TLAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:6
    plot(1:5,TRL(flyID,:),'Color','k','LineWidth',0.5);
end
for flyID = 1:5
    plot(1:5,TGL(flyID,:),'Color','k','LineWidth',0.6);    
end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,2)
hold on;
for vRng = 1:5
	scatter(zeros(6,1)+vRng,TRR(:,vRng),'k','filled');
    scatter(zeros(5,1)+vRng,TGR(:,vRng),'k','filled');
    TRAll = vertcat(TRR(:,vRng),TGR(:,vRng));
    line([vRng-1/2 vRng+1/2],[mean(TRAll(~isnan(TRAll))) mean(TRAll(~isnan(TRAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:6
    plot(1:5,TRR(flyID,:),'Color','k','LineWidth',0.5);
end
for flyID = 1:5
    plot(1:5,TGR(flyID,:),'Color','k','LineWidth',0.6);  
end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,3)
hold on;
for vRng = 1:5
	scatter(zeros(5,1)+vRng,WRL(:,vRng),'k');
    scatter(zeros(6,1)+vRng,WGL(:,vRng),'k');
    WLAll = vertcat(WRL(:,vRng),WGL(:,vRng));
    line([vRng-1/2 vRng+1/2],[mean(WLAll(~isnan(WLAll))) mean(WLAll(~isnan(WLAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:6
    plot(1:5,WGL(flyID,:),'Color','k','LineWidth',0.5);
end
for flyID = 1:5
    plot(1:5,WRL(flyID,:),'Color','k','LineWidth',0.6);  
end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,4)
hold on;
for vRng = 1:5
	scatter(zeros(5,1)+vRng,WRR(:,vRng),'k');
    scatter(zeros(6,1)+vRng,WGR(:,vRng),'k');
	WRAll = vertcat(WRR(:,vRng),WGR(:,vRng));
    line([vRng-1/2 vRng+1/2],[mean(WRAll(~isnan(WRAll))) mean(WRAll(~isnan(WRAll)))],'Color','k','LineWidth',2);
end
for flyID = 1:6
    plot(1:5,WGR(flyID,:),'Color','k','LineWidth',0.5);
end
for flyID = 1:5
    plot(1:5,WRR(flyID,:),'Color','k','LineWidth',0.6); 
end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

%%
load('D:\Imaging\2Color\20160903\DirDiff_60D05-VT8135.mat');
WGL = LdirDiffG;
WGR = RdirDiffG;
TRL = LdirDiffR;
TRR = RdirDiffR;
WRL = [];
WRR = [];
TGL = [];
TGR = [];

figure;
subplot(2,2,1)
hold on;
for vRng = 1:5
    scatter(zeros(3,1)+vRng,TRL(:,vRng),'k','filled');
%     scatter(zeros(5,1)+vRng,TGL(:,vRng),'k','filled');
    TLAll = TRL(:,vRng);
    line([vRng-1/2 vRng+1/2],[mean(TLAll(~isnan(TLAll))) mean(TLAll(~isnan(TLAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:3
    plot(1:5,TRL(flyID,:),'Color','k','LineWidth',0.5);
end
% for flyID = 1:5
%     plot(1:5,TGL(flyID,:),'Color','k','LineWidth',0.6);    
% end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,2)
hold on;
for vRng = 1:5
	scatter(zeros(3,1)+vRng,TRR(:,vRng),'k','filled');
%     scatter(zeros(5,1)+vRng,TGR(:,vRng),'k','filled');
    TRAll = TRR(:,vRng);
    line([vRng-1/2 vRng+1/2],[mean(TRAll(~isnan(TRAll))) mean(TRAll(~isnan(TRAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:3
    plot(1:5,TRR(flyID,:),'Color','k','LineWidth',0.5);
end
% for flyID = 1:5
%     plot(1:5,TGR(flyID,:),'Color','k','LineWidth',0.6);  
% end
ylim([-0.4 0.4]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,3)
hold on;
for vRng = 1:5
% 	scatter(zeros(5,1)+vRng,WRL(:,vRng),'k');
    scatter(zeros(3,1)+vRng,WGL(:,vRng),'k');
    WLAll = WGL(:,vRng);
    line([vRng-1/2 vRng+1/2],[mean(WLAll(~isnan(WLAll))) mean(WLAll(~isnan(WLAll)))],'Color','k','LineWidth',1);
end
for flyID = 1:3
    plot(1:5,WGL(flyID,:),'Color','k','LineWidth',0.5);
end
% for flyID = 1:5
%     plot(1:5,WRL(flyID,:),'Color','k','LineWidth',0.6);  
% end
ylim([-0.8 0.8]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

subplot(2,2,4)
hold on;
for vRng = 1:5
% 	scatter(zeros(5,1)+vRng,WRR(:,vRng),'k');
    scatter(zeros(3,1)+vRng,WGR(:,vRng),'k');
	WRAll = WGR(:,vRng);
    line([vRng-1/2 vRng+1/2],[mean(WRAll(~isnan(WRAll))) mean(WRAll(~isnan(WRAll)))],'Color','k','LineWidth',2);
end
for flyID = 1:3
    plot(1:5,WGR(flyID,:),'Color','k','LineWidth',0.5);
end
% for flyID = 1:5
%     plot(1:5,WRR(flyID,:),'Color','k','LineWidth',0.6); 
% end
ylim([-0.8 0.8]);
xlim([0.5 5.5]);
line([0 6],[0 0],'Color','k');

%% Look at fluorescence peak vs. rotational veloccity for 8135

% Specify the range of rotational velocities to consider
vRMin = 0;
vRMax = 720;
vRSpan = 30;

% Sort peaks by rotational velocity
[LPB.dark.dataR, LPB.dark.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(1:8),greenSpan(1:8), 'dark');
[RPB.dark.dataR, RPB.dark.dataG] = TurnSort(allFlyData, vRMin, vRMax, vRSpan, redSpan(9:16),greenSpan(9:16), 'dark');

% Look at the red channel
LRed = zeros(length(allFlyData),11);
RRed = zeros(length(allFlyData),11);

for flyID = 1:length(allFlyData)
    
    LRed(flyID,6) = mean(LPB.dark.dataR{flyID}.Stop);
    RRed(flyID,6) = mean(RPB.dark.dataR{flyID}.Stop);
    for vRng = 1:5
        
        LRed(flyID,6+vRng) = mean(LPB.dark.dataR{flyID}.CW{vRng});
        LRed(flyID,6-vRng) = mean(LPB.dark.dataR{flyID}.CCW{vRng});
        RRed(flyID,6+vRng) = mean(RPB.dark.dataR{flyID}.CW{vRng});
        RRed(flyID,6-vRng) = mean(RPB.dark.dataR{flyID}.CCW{vRng});
    end
    
end

figure;
subplot(2,2,1);
plot(LRed');
subplot(2,2,2);
plot(RRed');
subplot(2,2,3);
plot(LRed(:,7:11)'-fliplr(LRed(:,1:5))');
subplot(2,2,4);
plot(RRed(:,7:11)'-fliplr(RRed(:,1:5))');
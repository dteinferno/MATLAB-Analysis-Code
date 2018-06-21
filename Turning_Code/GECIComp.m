% Clear workspace
clear;
clc;
figure('units','normalized','outerposition',[0 0 0.5 1]);

% Set the model constants
ROIs = linspace(-pi, pi, 16);
tStep = 7/88; % imaging framerate

% Set the time constants.
% Taken from Fig. 7 - fig supp 3 from 
% Sensitive red protein calcium indicators for imaging neural activity
% Dana, H. et al., 2016 eLife
RtOn = 0.67/log(2);
RtOff = 0.52/log(2);
GtOn = 0.42/log(2);
GtOff = 0.33/log(2);

for shft = 1:3
bumpShift = (shft-1)*pi/32;

for velBump = 1:8
    
bumpSpeed = 10*(velBump) * pi/180; % speed of the bump around the ring
tAll = [0:tStep:10*ceil(2*pi/bumpSpeed)]; % time


% Convolute the rising and falling time constants with the Ca signal
% (assumed to be a cos shape around the ring) and plot
for ROI=1:16
    RbumpShape = 0.5*cos(ROIs(ROI)-bumpSpeed*tAll)+0.5;
    RbumpShapeUp = zeros(1,length(RbumpShape));
    RgoingUp = find(diff(RbumpShape)>0);
    RbumpShapeUp(RgoingUp) = RbumpShape(RgoingUp);
    RbumpShapeDown = zeros(length(RbumpShape),1);
    RgoingDown = find(diff(RbumpShape)<0);
    RbumpShapeDown(RgoingDown) = RbumpShape(RgoingDown);

    RbumpUp = conv(RbumpShapeUp,exp(-tAll/RtOn));
    RbumpDown = conv(RbumpShapeDown,exp(-tAll/RtOff));
    Rbump(ROI,:) = RbumpUp + RbumpDown;
    
    GbumpShape = 0.5*cos(ROIs(ROI)-bumpShift-bumpSpeed*tAll)+0.5;
    GbumpShapeUp = zeros(1,length(GbumpShape));
    GgoingUp = find(diff(GbumpShape)>0);
    GbumpShapeUp(GgoingUp) = GbumpShape(GgoingUp);
    GbumpShapeDown = zeros(length(GbumpShape),1);
    GgoingDown = find(diff(GbumpShape)<0);
    GbumpShapeDown(GgoingDown) = GbumpShape(GgoingDown);
    
    GbumpUp = conv(GbumpShapeUp,exp(-tAll/GtOn));
    GbumpDown = conv(GbumpShapeDown,exp(-tAll/GtOff));
    Gbump(ROI,:) = GbumpUp + GbumpDown;
      
end

if velBump == 3
    subplot(3,3,1+(shft-1)*3);
    hold on;
    title(strcat('v_r = ',num2str(velBump*10)));
    plot(ROIs,0.5*cos(ROIs-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
    plot(ROIs,0.5*cos(ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
    plot(ROIs,Rbump(:,round(length(tAll/2)))./max(Rbump(:,round(length(tAll/2)))),'r','LineWidth',2);
    plot(ROIs,Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2)))),'g','LineWidth',2);
end

RturnsR = Rbump(:,round(length(tAll/2)))./max(Rbump(:,round(length(tAll/2))));
RturnsG = Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2))));
angsraw = (1:16)*2*pi/16+pi;
angsraw = angsraw';
PVAdiff1(velBump) = circ_mean(angsraw,RturnsG)-circ_mean(angsraw,RturnsR);
PVAdiff1(velBump) = min(PVAdiff1(velBump),2*pi-PVAdiff1(velBump));
       
clear Gbump RturnsR RturnsG;

bumpShift = -bumpShift;
for ROI=1:16
    GbumpShape = 0.5*cos(ROIs(ROI)-bumpShift-bumpSpeed*tAll)+0.5;
    GbumpShapeUp = zeros(1,length(GbumpShape));
    GgoingUp = find(diff(GbumpShape)>0);
    GbumpShapeUp(GgoingUp) = GbumpShape(GgoingUp);
    GbumpShapeDown = zeros(length(GbumpShape),1);
    GgoingDown = find(diff(GbumpShape)<0);
    GbumpShapeDown(GgoingDown) = GbumpShape(GgoingDown);
    
    GbumpUp = conv(GbumpShapeUp,exp(-tAll/GtOn));
    GbumpDown = conv(GbumpShapeDown,exp(-tAll/GtOff));
    Gbump(ROI,:) = GbumpUp + GbumpDown;
end

if velBump==3
    plot(ROIs,0.5*cos(ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
    plot(ROIs,Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2)))),'g','LineWidth',2);

    set(gca,'FontSize',16);
    xlabel('EB location (rad)');
    xlim([-pi pi]);
end

RturnsR = Rbump(:,round(length(tAll/2)))./max(Rbump(:,round(length(tAll/2))));
RturnsG = Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2))));
PVAdiff2(velBump) = circ_mean(angsraw,RturnsR)-circ_mean(angsraw,RturnsG);
PVAdiff2(velBump) = min(PVAdiff2(velBump),2*pi-PVAdiff2(velBump));

clear Gbump Rbump;
end

for ln = 1:length(PVAdiff1)/2
    subplot(3,3,2+(shft-1)*3);
    hold on;
    if ln == 1
        rectangle('Position',[-0.15 min(0,-180/pi*PVAdiff1(2)) 0.3  abs(180/pi*PVAdiff1(2))]);
    else
        rectangle('Position',[ln-1-0.15 -180/pi*max(PVAdiff1(2*ln),PVAdiff1(2*ln-2)) 0.3 180/pi*abs(PVAdiff1(2*ln)-PVAdiff1(2*ln-2))]);
    end
    
    subplot(3,3,3+(shft-1)*3);
    hold on;
    if ln == 1
        rectangle('Position',[-0.15 min(0,180/pi*PVAdiff2(2)) 0.3  abs(180/pi*PVAdiff2(2))]);
    else
        rectangle('Position',[ln-1-0.15 180/pi*PVAdiff2(2*ln-2) 0.3 180/pi*abs(PVAdiff2(2*ln)-PVAdiff2(2*ln-2))]);
    end
end

clear PVAdiff1 PVAdiff2;

subplot(3,3,2+(shft-1)*3);
load('FlipPVADiff');
flipPVAdiff = -PVAdiff;
scatter(zeros(5,1),180/pi*flipPVAdiff(:,1),10,'filled','k');
line([-0.1 0.1],[180/pi*mean(flipPVAdiff(:,1)) 180/pi*mean(flipPVAdiff(:,1))],'Color','k','LineWidth',1);
scatter(zeros(5,1)+1,180/pi*flipPVAdiff(:,2),10,'filled','k');
line([0.9 1.1],[180/pi*mean(flipPVAdiff(:,2)) 180/pi*mean(flipPVAdiff(:,2))],'Color','k','LineWidth',1);
scatter(zeros(5,1)+2,180/pi*flipPVAdiff(:,3),10,'filled','k');
line([1.9 2.1],[180/pi*mean(flipPVAdiff(:,3)) 180/pi*mean(flipPVAdiff(:,3))],'Color','k','LineWidth',1);
scatter(zeros(5,1)+3,180/pi*flipPVAdiff(:,4),10,'filled','k');
line([2.9 3.1],[180/pi*mean(flipPVAdiff(:,4)) 180/pi*mean(flipPVAdiff(:,4))],'Color','k','LineWidth',1);

set(gca,'FontSize',14);
ylabel('PVA difference (deg)');
set(gca,'XTickLabel',{'0-20','20-40','40-60','60-80'});
xlabel('v_r bins (deg/sec)');
xlim([-0.5 3.5]);
ylim([-20 50]);

subplot(3,3,3+(shft-1)*3);
load('PVAdiff');
scatter(zeros(10,1),180/pi*PVAdiff(:,1),10,'filled','k');
line([-0.1 0.1],[180/pi*mean(PVAdiff(:,1)) 180/pi*mean(PVAdiff(:,1))],'Color','k','LineWidth',1);
scatter(zeros(10,1)+1,180/pi*PVAdiff(:,2),10,'filled','k');
line([0.9 1.1],[180/pi*mean(PVAdiff(:,2)) 180/pi*mean(PVAdiff(:,2))],'Color','k','LineWidth',1);
scatter(zeros(10,1)+2,180/pi*PVAdiff(:,3),10,'filled','k');
line([1.9 2.1],[180/pi*mean(PVAdiff(:,3)) 180/pi*mean(PVAdiff(:,3))],'Color','k','LineWidth',1);
scatter(zeros(10,1)+3,180/pi*PVAdiff(:,4),10,'filled','k');
line([2.9 3.1],[180/pi*mean(PVAdiff(:,4)) 180/pi*mean(PVAdiff(:,4))],'Color','k','LineWidth',1);

set(gca,'FontSize',14);
ylabel('PVA difference (deg)');
set(gca,'XTickLabel',{'0-20','20-40','40-60','60-80'});
xlabel('v_r bins (deg/sec)');
xlim([-0.5 3.5]);
ylim([-20 50]);

end
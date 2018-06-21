%% Given bump offset vs vRot function and time constants,
% compare calculated and measured delta PVA vs. vRot

% Clear workspace
clear;
clc;

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

% general graphics, this will apply to any figure you open (groot is the default figure object).
% I have this in my startup.m file, so I don't have to retype these things whenever plotting a new fig.
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 8, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);

GECIComp = figure('units','normalized','outerposition',[0 0 0.5 1]);

% Set the model constants
ROIs = linspace(-pi, pi, 16);
FWHM = pi/2;
tStep = 0.0876; % imaging framerate

invert = 3;
if invert == 1
    % Set the time constants.
    % Taken from Fig. 7 - fig supp 3 from 
    % Sensitive red protein calcium indicators for imaging neural activity
    % Dana, H. et al., 2016 eLife
    RtOn = 0.67/log(2);
    RtOff = 0.52/log(2);
    GtOn = 0.42/log(2);
    GtOff = 0.33/log(2);    
elseif invert == 2
    RtOn = 0.2/log(2);
    RtOff = 0.45/log(2);
    GtOn = 0.16/log(2);
    GtOff = 0.42/log(2);
else
    % Alternative time constants taken from Fig 2 for 10 Hz
    RtOn = 0.078/log(2);
    RtOff = 0.43/log(2);
    % Time constants taken from Ultrasensitive fluorescent proteins for
    % imaging neuronal activity for 10 Hz
    GtOn = 0.092/log(2);
    GtOff = 0.44/log(2);
end

% Specify the velocities of interest
vRs = [0:15:180];

% Specify the bump shift for each vR value
bumpShiftAll = zeros(4,length(vRs));
bumpShiftAll(1,:) = 0;
subplot(5,3,4);
plot(vRs,180/pi*bumpShiftAll(1,:),'k');
xlim([0 180]);
ylim([-20 45]);
xlabel('vrot (o/s)')
ylabel('bump offset (o)')
set(gca,'FontSize',10);
set(gca,'XTick',[0:30:180]);
box off;

bumpShiftAll(2,:) = pi/16;
subplot(5,3,7);
plot(vRs,180/pi*bumpShiftAll(2,:),'k');
xlim([0 180]);
ylim([-20 20]);
xlabel('vrot (o/s)')
ylabel('bump offset (o)')
set(gca,'FontSize',10);
set(gca,'XTick',[0:30:180]);
box off;

bumpShiftAll(3,:) = vRs* (3*pi/32/165);
subplot(5,3,10);
plot(vRs,180/pi*bumpShiftAll(3,:),'k');
xlim([0 180]);
ylim([-20 45]);
xlabel('vrot (o/s)')
ylabel('bump offset (o)')
set(gca,'FontSize',10);
set(gca,'XTick',[0:30:180]);
box off;

bumpShiftAll(4,:) = 3*pi/32*(1-exp(-vRs/45));
subplot(5,3,13);
plot(vRs,180/pi*bumpShiftAll(4,:),'k');
xlim([0 180]);
ylim([-20 45]);
xlabel('vrot (o/s)')
ylabel('bump offset (o)')
set(gca,'FontSize',10);
set(gca,'XTick',[0:30:180]);
box off;


% for i=1:length(vRs)
%     bumpShiftAll(3,i) = pi/16*1/(1+exp(-0.1*(vRs(i)-90)));
% end

% Step through the different shift models
for shft = 1:4

    PVAdiff1(1) = bumpShiftAll(shft,1);
    PVAdiff2(1) = bumpShiftAll(shft,1);
    % Advance through the velocities
    for velBump = 2:length(vRs)

        bumpShift = bumpShiftAll(shft,velBump); % bump shift
        bumpSpeed = vRs(velBump) * pi/180; % speed of the bump around the ring
        if bumpSpeed == 0
            tAll = [0:tStep:10];
        else
            tAll = [0:tStep:10*ceil(2*pi/bumpSpeed)]; % Look at the bumps after 10 rotations
        end


        % Convolute the rising and falling time constants with the Ca signal
        % (assumed to be a cos shape around the ring) and plot
        for ROI=1:16
%             RbumpShape = 0.5*cos(ROIs(ROI)-bumpSpeed*tAll)+0.5;
            RbumpShape = exp(-(mod((ROIs(ROI)-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Rbump(ROI,:) = conv(RbumpShape,abs(exp(-tAll/RtOn)-exp(-tAll/RtOff)));
            
%             GbumpShape = 0.5*cos(ROIs(ROI)-bumpShift-bumpSpeed*tAll)+0.5;
            GbumpShape = exp(-(mod((ROIs(ROI)-bumpShift-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Gbump(ROI,:) = conv(GbumpShape,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));

        end

        if ( velBump == 2 || velBump == 7 || velBump == 12 ) & shft == 3
            subplot(5,3,(velBump-2)/5+1);
            hold on;
            title(strcat('vrot = ',num2str(velBump*15-15),'o/s'));
            RBPlot = Rbump(:,round(length(tAll/2)))./max(Rbump(:,round(length(tAll/2))));
            RBMax = find(RBPlot == max(RBPlot));
            RBPlot(1:end-1) = circshift(RBPlot(1:end-1),8-RBMax);
            RBPlot(end) = RBPlot(1);
            GBPlot = Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2))));
            GBPlot(1:end-1) = circshift(GBPlot(1:end-1),8-RBMax);
            GBPlot(end) = GBPlot(1);
            %             plot(ROIs,0.5*cos(ROIs-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
            plot(ROIs,circshift(exp(-(mod((ROIs-bumpSpeed*tAll(round(length(tAll/2)))),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2))',8-RBMax),'k');
%             plot(ROIs,0.5*cos(ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
            plot(ROIs,circshift(exp(-(mod((ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2)))),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2))',8-RBMax),'k');


            plot(ROIs,RBPlot,'m','LineWidth',2);
            plot(ROIs,GBPlot,'g','LineWidth',2);
            xlim([-pi pi]);
            set(gca,'FontSize',10);
            xlabel('EB position (rad)');
            ylabel('bump amplitude (A.U.)');
        end

        RturnsR = Rbump(:,floor(length(tAll/2)))./max(Rbump(:,floor(length(tAll/2))));
        RturnsG = Gbump(:,floor(length(tAll/2)))./max(Gbump(:,floor(length(tAll/2))));
        angsraw = (1:16)*2*pi/16+pi;
        angsraw = angsraw';
        PVAdiff1(velBump) = circ_mean(angsraw,RturnsG)-circ_mean(angsraw,RturnsR);
        if PVAdiff1(velBump) > pi
            PVAdiff1(velBump) = PVAdiff1(velBump)-2*pi;
        elseif PVAdiff1(velBump) < -pi
            PVAdiff1(velBump) = PVAdiff1(velBump)+2*pi;
        end

        clear Gbump RturnsR RturnsG;
        
        for ROI=1:16
%             GbumpShape = 0.5*cos(ROIs(ROI)-bumpShift-bumpSpeed*tAll)+0.5;
            GbumpShape = exp(-(mod((ROIs(ROI)+bumpShift-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Gbump(ROI,:) = conv(GbumpShape,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        end

%         if velBump == 5
% %             plot(ROIs,0.5*cos(ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
% %             plot(ROIs,exp(-(mod((ROIs+bumpShift-bumpSpeed*tAll(round(length(tAll/2)))),2*pi)-pi).^2/(2*(pi/2)^2)),'g');
% %             plot(ROIs,Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2)))),':g','LineWidth',2);
% 
%             set(gca,'FontSize',16);
%             xlabel('EB location (rad)');
%             xlim([-pi pi]);
%         end

        RturnsR = Rbump(:,floor(length(tAll/2)))./max(Rbump(:,floor(length(tAll/2))));
        RturnsG = Gbump(:,floor(length(tAll/2)))./max(Gbump(:,floor(length(tAll/2))));
        PVAdiff2(velBump) = circ_mean(angsraw,RturnsR)-circ_mean(+angsraw,RturnsG);
        if PVAdiff2(velBump) > pi
            PVAdiff2(velBump) = PVAdiff2(velBump)-2*pi;
        elseif PVAdiff2(velBump) < -pi
            PVAdiff2(velBump) = PVAdiff2(velBump)+2*pi;
        end
        
%         if ( velBump == 2 || velBump == 7 || velBump == 12 ) & shft == 2
%             subplot(5,3,(velBump-2)/5+1);
%             hold on;
%             title(strcat('vrot = ',num2str(velBump*15-15),'o/s'));
%             RBPlot = Rbump(:,round(length(tAll/2)))./max(Rbump(:,round(length(tAll/2))));
%             RBMax = find(RBPlot == max(RBPlot));
%             RBPlot(1:end-1) = circshift(RBPlot(1:end-1),8-RBMax);
%             RBPlot(end) = RBPlot(1);
%             GBPlot = Gbump(:,round(length(tAll/2)))./max(Gbump(:,round(length(tAll/2))));
%             GBPlot(1:end-1) = circshift(GBPlot(1:end-1),8-RBMax);
%             GBPlot(end) = GBPlot(1);
%             %             plot(ROIs,0.5*cos(ROIs-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
%             plot(ROIs,circshift(exp(-(mod((ROIs-bumpSpeed*tAll(round(length(tAll/2)))),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2))',8-RBMax),'k');
% %             plot(ROIs,0.5*cos(ROIs-bumpShift-bumpSpeed*tAll(round(length(tAll/2))))+0.5,'k');
%             plot(ROIs,circshift(exp(-(mod((ROIs+bumpShift-bumpSpeed*tAll(round(length(tAll/2)))),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2))',8-RBMax),'k');
% 
% 
%             plot(ROIs,RBPlot,'m','LineWidth',2);
%             plot(ROIs,GBPlot,'g','LineWidth',2);
%             xlim([-pi pi]);
%             set(gca,'FontSize',10);
%             xlabel('EB position (rad)');
%             ylabel('bump amplitude (A.U.)');
%         end

        clear Gbump Rbump;
    end

    PVAdiff1 = 180/pi*PVAdiff1;
    PVAdiff2 = 180/pi*PVAdiff2;    
    for ln = 1:floor(length(PVAdiff1)/2)
        subplot(5,3,5+(shft-1)*3);
        hold on;
        if ln == 1
            rectangle('Position',[-0.15 min(PVAdiff1(2*ln-1),PVAdiff1(2*ln+1)) 0.3  abs(PVAdiff1(3)-PVAdiff1(1))],'FaceColor',[0.8 0.8 0.8]);
        else
            rectangle('Position',[ln-1-0.15 min(PVAdiff1(2*ln-1),PVAdiff1(2*ln+1)) 0.3 abs(PVAdiff1(2*ln+1)-PVAdiff1(2*ln-1))],'FaceColor',[0.8 0.8 0.8]);
        end

        subplot(5,3,6+(shft-1)*3);
        hold on;
        if ln == 1
            rectangle('Position',[-0.15 min(PVAdiff2(2*ln-1),PVAdiff2(2*ln+1)) 0.3  abs(PVAdiff2(3)-PVAdiff2(1))],'FaceColor',[0.8 0.8 0.8]);
        else
            rectangle('Position',[ln-1-0.15 min(PVAdiff2(2*ln-1),PVAdiff2(2*ln+1)) 0.3 abs(PVAdiff2(2*ln+1)-PVAdiff2(2*ln-1))],'FaceColor',[0.8 0.8 0.8]);
        end
    end

    clear PVAdiff1 PVAdiff2;

    subplot(5,3,5+(shft-1)*3);
    load('D:\Imaging\2Color\Stats\EB\PVADiffR60D05G37F06');
    PVAdiff = adjPVAdiff;
    for vPlt = 1:6
        scatter(zeros(10,1)+vPlt-1,PVAdiff(:,vPlt),10,'filled','k');
%         line([-0.1 0.1],[mean(PVAdiff(:,vPlt)) mean(PVAdiff(:,vPlt))],'Color','k','LineWidth',1);
    end
    set(gca,'FontSize',10);
    ylabel('PVA difference (o)');
    set(gca,'XTick',[0:6]);
    set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
    xlabel('vrot bins (o/s)');
    xlim([-0.5 5.5]);
    ylim([-20 45]);

    subplot(5,3,6+(shft-1)*3);
    load('D:\Imaging\2Color\Stats\EB\PVADiffR37F06G60D05')
    flipPVAdiff = adjPVAdiff;
    for vPlt = 1:6
        scatter(zeros(5,1)+vPlt-1,flipPVAdiff(:,vPlt),10,'filled','k');
%         line([-0.1 0.1],[mean(flipPVAdiff(:,vPlt)) mean(flipPVAdiff(:,vPlt))],'Color','k','LineWidth',1);
    end
    set(gca,'FontSize',10);
    ylabel('PVA difference (o)');
    set(gca,'XTick',[0:6]);
    set(gca,'XTickLabel',{'0-30','30-60','60-90','90-120','120-150','150-180'});
    xlabel('vrot bins (o/s)');
    xlim([-0.5 5.5]);
    ylim([-20 45]);
    
end

set(GECIComp,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
if invert == 1
    print(GECIComp,'D:\Imaging\2Color\GECIComp','-dpdf');
elseif invert == 2
    print(GECIComp,'D:\Imaging\2Color\GECIComp_Fit','-dpdf');
else
    print(GECIComp,'D:\Imaging\2Color\GECIComp_Vert','-dpdf');
end
%% Given bump offset vs vRot function and time constants,
% compare calculated and measured delta PVA vs. vRot in the PB

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
ROIs = linspace(-pi, pi, 8);
FWHM = pi/2;
tStep = 0.0876; % imaging framerate

RtOnAll = [0.67 0.2 0.078]./log(2);
RtOffAll = [0.52 0.45 0.43]./log(2);
GtOnAll = [0.42 0.16 0.092]./log(2);
GtOffAll = [0.33 0.42 0.44]./log(2);

% Specify the velocities of interest
vRs = [0:15:150];

% Specify the bump shift for each vR value
bumpShiftAll = vRs*1/165*15*(8/360);

% Load the data
load('D:\Imaging\2Color\20160621\PVADiff_60D05-37F062.mat')
LWmT = LPVAdiff;
RWmT = RPVAdiff;
load('D:\Imaging\2Color\20160622\PVADiff_37F06-60D052.mat')
LTmW= LPVAdiff;
RTmW = RPVAdiff;

% Set the data plotting parameters
vRMin = 0;  
vRMax = 720;
vRSpan = 30;
plotRange = 150/vRSpan;

% Step through the different shift models
for invert = 1:3

    RtOn = RtOnAll(invert);
    RtOff = RtOffAll(invert);
    GtOn = GtOnAll(invert);
    GtOff = GtOffAll(invert);
    
    PVAdiff1(1) = bumpShiftAll(1);
    PVAdiff2(1) = bumpShiftAll(1);
    % Advance through the velocities
    for velBump = 2:length(vRs)

        bumpShift = bumpShiftAll(velBump); % bump shift
        bumpSpeed = vRs(velBump) * pi/180; % speed of the bump around the ring
        if bumpSpeed == 0
            tAll = [0:tStep:10];
        else
            tAll = [0:tStep:10*ceil(2*pi/bumpSpeed)]; % Look at the bumps after 10 rotations
        end


        % Convolute the rising and falling time constants with the Ca signal
        % (assumed to be a cos shape around the ring) and plot
        for ROI=1:8
            RbumpShape = exp(-(mod((ROIs(ROI)-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Rbump(ROI,:) = conv(RbumpShape,abs(exp(-tAll/RtOn)-exp(-tAll/RtOff)));
            
            GbumpShape = exp(-(mod((ROIs(ROI)-bumpShift-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Gbump(ROI,:) = conv(GbumpShape,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));

        end

        RturnsR = Rbump(:,floor(length(tAll/2)))./max(Rbump(:,floor(length(tAll/2))));
        RturnsG = Gbump(:,floor(length(tAll/2)))./max(Gbump(:,floor(length(tAll/2))));
        angsraw = (1:8)*2*pi/8+pi;
        angsraw = angsraw';
        PVAdiff1(velBump) = circ_mean(angsraw,RturnsG)-circ_mean(angsraw,RturnsR);
        if PVAdiff1(velBump) > pi
            PVAdiff1(velBump) = PVAdiff1(velBump)-2*pi;
        elseif PVAdiff1(velBump) < -pi
            PVAdiff1(velBump) = PVAdiff1(velBump)+2*pi;
        end

        clear Gbump RturnsR RturnsG;
        
        for ROI=1:8
            GbumpShape = exp(-(mod((ROIs(ROI)+bumpShift-bumpSpeed*tAll),2*pi)-pi).^2/(2*(FWHM/(2*sqrt(2*log(2))))^2));
            Gbump(ROI,:) = conv(GbumpShape,abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
        end

        RturnsR = Rbump(:,floor(length(tAll/2)))./max(Rbump(:,floor(length(tAll/2))));
        RturnsG = Gbump(:,floor(length(tAll/2)))./max(Gbump(:,floor(length(tAll/2))));
        PVAdiff2(velBump) = circ_mean(angsraw,RturnsR)-circ_mean(+angsraw,RturnsG);
        if PVAdiff2(velBump) > pi
            PVAdiff2(velBump) = PVAdiff2(velBump)-2*pi;
        elseif PVAdiff2(velBump) < -pi
            PVAdiff2(velBump) = PVAdiff2(velBump)+2*pi;
        end

        clear Gbump Rbump;
    end

    PVAdiff1 = 4/pi*PVAdiff1;
    PVAdiff2 = 4/pi*PVAdiff2;    
    for ln = 1:floor(length(PVAdiff1)/2)
        subplot(3,2,invert*2-1);
        hold on;
        if ln == 1
            rectangle('Position',[vRs(1) min(PVAdiff1(2*ln-1),PVAdiff1(2*ln+1)) 30  abs(PVAdiff1(3)-PVAdiff1(1))],'FaceColor',[0.8 0.8 0.8]);
        else
            rectangle('Position',[vRs(2*ln-1) min(PVAdiff1(2*ln-1),PVAdiff1(2*ln+1)) 30 abs(PVAdiff1(2*ln+1)-PVAdiff1(2*ln-1))],'FaceColor',[0.8 0.8 0.8]);
        end

        subplot(3,2,2*invert);
        hold on;
        if ln == 1
            rectangle('Position',[vRs(1) min(PVAdiff2(2*ln-1),PVAdiff2(2*ln+1)) 30  abs(PVAdiff2(3)-PVAdiff2(1))],'FaceColor',[0.8 0.8 0.8]);
        else
            rectangle('Position',[vRs(2*ln-1) min(PVAdiff2(2*ln-1),PVAdiff2(2*ln+1)) 30 abs(PVAdiff2(2*ln+1)-PVAdiff2(2*ln-1))],'FaceColor',[0.8 0.8 0.8]);
        end
    end

    clear PVAdiff1 PVAdiff2;
    
    subplot(3,2,2*invert-1);
    hold on;
    for bnz = 1:(length(LWmT)-1)/2
        LNow = -LTmW(:,bnz);
        LNow(find(LNow==0)) = [];
        scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'k','filled');
        RNow = -RTmW(:,bnz);
        RNow(find(RNow==0)) = [];
        scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'k','filled');
    end
    for bnz = (length(LWmT)+1)/2+1:length(LWmT)
        LNow = LTmW(:,bnz);
        LNow(find(LNow==0)) = [];
        scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'k','filled');
        RNow = RTmW(:,bnz);
        RNow(find(RNow==0)) = [];
        scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'k','filled');
    end
    line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
    xlim([-3*vRSpan/4-vRSpan/10 150]);
    ylim([-1.5 2.5]);
    ylabel('offset (# of glom.)');

    subplot(3,2,2*invert);
    hold on;
    for bnz = 1:(length(LWmT)-1)/2
        LNow = LWmT(:,bnz);
        LNow(find(LNow==0)) = [];
        scatter(zeros(length(LNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2-vRSpan/10,LNow,'k','filled');
        RNow = RWmT(:,bnz);
        RNow(find(RNow==0)) = [];
        scatter(zeros(length(RNow),1)-(bnz-1)*vRSpan+plotRange*vRSpan-vRSpan/2+vRSpan/10,RNow,'k','filled');
    end
    for bnz = (length(LWmT)+1)/2+1:length(LWmT)
        LNow = -LWmT(:,bnz);
        LNow(find(LNow==0)) = [];
        scatter(zeros(length(LNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2+vRSpan/10,LNow,'k','filled');
        RNow = -RWmT(:,bnz);
        RNow(find(RNow==0)) = [];
        scatter(zeros(length(RNow),1)+(bnz-1)*vRSpan-plotRange*vRSpan-vRSpan/2-vRSpan/10,RNow,'k','filled');
    end
    line([-3*vRSpan/4-vRSpan/10 vRSpan*plotRange],[0 0],'Color','k');
    xlim([-3*vRSpan/4-vRSpan/10 150]);
    ylim([-1.5 2.5]); 
end

subplot(3,2,1);
title('P-ENs: green, E-PGs: red');
subplot(3,2,2);
title('P-ENs: red, E-PGs: green');
subplot(3,2,5);
xlabel('vRot (o/s)');
subplot(3,2,6);
xlabel('vRot (o/s)');

set(GECIComp,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(GECIComp,'D:\Imaging\2Color\GECIComp_PB2','-dpdf');
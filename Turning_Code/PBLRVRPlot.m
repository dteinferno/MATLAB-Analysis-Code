flyID = 3;

figure;

for type = 1:4
    if type == 1
        data = LPB.dark.dataR{flyID};
    elseif type == 2
        data = RPB.dark.dataR{flyID};
    elseif type == 3
        data = LPB.dark.dataG{flyID};
    elseif type == 4
        data = RPB.dark.dataG{flyID};
    end
        
    % Specify the range of rotational velocities to consider
    vRMin = 0*pi/180;
    vRMax = 720*pi/180;
    vRSpan = 30;
    plotMax = 150/vRSpan;

    CWStdLow = zeros(1,plotMax);
    CWMean = zeros(1,plotMax);
    CWStdHigh = zeros(1,plotMax);
    CCWStdLow = zeros(1,plotMax);
    CCWMean = zeros(1,plotMax);
    CCWStdHigh = zeros(1,plotMax);

    for vRng = 1:plotMax
        CWStdLow(vRng) = mean(data.CW{vRng})-std(data.CW{vRng});
        CWMean(vRng) = mean(data.CW{vRng});
        CWStdHigh(vRng) = mean(data.CW{vRng})+std(data.CW{vRng});

        CCWStdLow(vRng) = mean(data.CCW{vRng})-std(data.CCW{vRng});
        CCWMean(vRng) = mean(data.CCW{vRng});
        CCWStdHigh(vRng) = mean(data.CCW{vRng})+std(data.CCW{vRng});
    end

    StopStdLow = mean(data.Stop)-std(data.Stop);
    StopMean = mean(data.Stop);
    StopStdHigh = mean(data.Stop)+std(data.Stop);

    subplot(2,2,type);
    hold on;
    scatter(0,StopMean,'k');
    line([0 0],[StopStdLow StopStdHigh],'Color','k');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CWStdLow,'LineWidth',1,'Color','c');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CWMean,'LineWidth',2,'Color','c');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CWStdHigh,'LineWidth',1,'Color','c');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CCWStdLow,'LineWidth',1,'Color','m');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CCWMean,'LineWidth',2,'Color','m');
    plot(-vRSpan/2+vRSpan*[1:plotMax],CCWStdHigh,'LineWidth',1,'Color','m');

end

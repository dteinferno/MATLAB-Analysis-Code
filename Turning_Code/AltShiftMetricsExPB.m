figure;
RPlt = LPBDataRAllCCW;
GPlt = LPBDataGAllCCW;



subplot(2,2,1);
hold on;
RPlt = mean(LPBDataRAllCCW,2);
RPltSD = std(LPBDataRAllCCW,[],2);

RPlt = (RPlt-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
RPltSD = (RPltSD-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
plot(1:9, RPlt,'r');
plot(1:9, RPlt+RPltSD,'r');
plot(1:9, RPlt-RPltSD,'r');
line([4+circ_mean(angsraw,RPlt(1:8))*num_ROIs/(2*pi) 4+circ_mean(angsraw,RPlt(1:8))*num_ROIs/(2*pi)],...
    [0 1],'color','r');

GPlt = mean(LPBDataGAllCCW,2);
GPltSD = std(LPBDataGAllCCW,[],2);
GPlt = (GPlt-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
GPltSD = (GPltSD-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
plot(1:9, GPlt,'g');
plot(1:9, GPlt+GPltSD,'g');
plot(1:9, GPlt-GPltSD,'g');
line([5+circ_mean(angsraw,GPlt(2:9))*num_ROIs/(2*pi) 5+circ_mean(angsraw,GPlt(2:9))*num_ROIs/(2*pi)],...
    [0 1],'color','g');
ylim([-0.5 1.5]);

subplot(2,2,2);
hold on;
for i = 1:5
    RPlt = LPBDataRAllCCW(:,i);
    RPlt = (RPlt-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
    plot(1:9, RPlt,'r');
    line([4+circ_mean(angsraw,RPlt(1:8))*num_ROIs/(2*pi) 4+circ_mean(angsraw,RPlt(1:8))*num_ROIs/(2*pi)],...
    [0 1],'color','r');

    GPlt = LPBDataGAllCCW(:,i);
    GPlt = (GPlt-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
    plot(1:9, GPlt,'g');
    line([5+circ_mean(angsraw,GPlt(2:9))*num_ROIs/(2*pi) 5+circ_mean(angsraw,GPlt(2:9))*num_ROIs/(2*pi)],...
    [0 1],'color','g');
end
ylim([-0.5 1.5]);

subplot(2,2,3);
hold on;
RPlt = mean(LPBDataRAllCCW,2);
RPltSD = std(LPBDataRAllCCW,[],2);

RPlt = (RPlt-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
RPltSD = (RPltSD-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
plot(1:9, RPlt,'r');
plot(1:9, RPlt+RPltSD,'r');
plot(1:9, RPlt-RPltSD,'r');
line([find(RPlt == max(RPlt)) find(RPlt == max(RPlt))],...
    [0 1],'color','r');

GPlt = mean(LPBDataGAllCCW,2);
GPltSD = std(LPBDataGAllCCW,[],2);
GPlt = (GPlt-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
GPltSD = (GPltSD-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
plot(1:9, GPlt,'g');
plot(1:9, GPlt+GPltSD,'g');
plot(1:9, GPlt-GPltSD,'g');
line([find(GPlt == max(GPlt)) find(GPlt == max(GPlt))],...
    [0 1],'color','g');
ylim([-0.5 1.5]);

subplot(2,2,4);
hold on;
for i = 1:5
    RPlt = LPBDataRAllCCW(:,i);
    RPlt = (RPlt-min(RPlt(1:8)))./(max(RPlt(1:8)-min(RPlt(1:8))));
    plot(1:9, RPlt,'r');
    line([find(RPlt == max(RPlt)) find(RPlt == max(RPlt))],...
    [0 1],'color','r');

    GPlt = LPBDataGAllCCW(:,i);
    GPlt = (GPlt-min(GPlt(2:9)))./(max(GPlt(2:9)-min(GPlt(2:9))));
    plot(1:9, GPlt,'g');
    line([find(GPlt == max(GPlt)) find(GPlt == max(GPlt))],...
    [0 1],'color','g');
end
ylim([-0.5 1.5]);
function TurnPlot(data, Min, Max, Span, plotMax)
%function to make scatter plots for clockwise and counterclockwise turns

%input:
%data = the data to be plotted
%Min = the minimum velocity bin (deg)
%Max = the maximum velocity bin (deg)
%Span = the span of each bin (deg)

%output:

% Specify the range of rotational velocities to consider
vRMin = Min*pi/180;
vRMax = Max*pi/180;
vRSpan = Span*pi/180;

for vRng = 1:plotMax
    if ~isempty(data.CW{vRng});
%         scatter(vRSpan/2+vRSpan*(vRng-1),median(data.CW{vRng}),'c');
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CW{vRng}),'c','filled');
%         line([vRSpan/2+vRSpan*(vRng-1) vRSpan/2+vRSpan*(vRng-1)],...
%             [quantile(data.CW{vRng},0.25) ...
%             quantile(data.CW{vRng},0.75)],...
%             'Color','c','LineWidth',1);
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CW{vRng})-std(data.CW{vRng}),10,'c');
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CW{vRng})+std(data.CW{vRng}),10,'c');

    end
    if ~isempty(data.CCW{vRng});
%         scatter(vRSpan/2+vRSpan*(vRng-1),median(data.CCW{vRng}),'m');
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CCW{vRng}),'m','filled');
%         line([vRSpan/2+vRSpan*(vRng-1) vRSpan/2+vRSpan*(vRng-1)],...
%             [quantile(data.CCW{vRng},0.25) ...
%             quantile(data.CCW{vRng},0.75)],...
%             'Color','m','LineWidth',1);
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CCW{vRng})-std(data.CCW{vRng}),10,'m');
        scatter(vRSpan/2+vRSpan*(vRng-1),mean(data.CCW{vRng})+std(data.CCW{vRng}),10,'m');
    end
end
% scatter(0,median(data.Stop),'k');
scatter(0,mean(data.Stop),'k');
% line([0 0],...
%             [quantile(data.Stop,0.25) ...
%             quantile(data.Stop,0.75)],...
%             'Color','k','LineWidth',1);
scatter(0,mean(data.Stop)-std(data.Stop),10,'k');
scatter(0,mean(data.Stop)+std(data.Stop),10,'k');
    
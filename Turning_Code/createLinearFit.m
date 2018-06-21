
function [fitresult, gof] = createLinearFit(x, y, w)
 
%  Input:
%       x: angular velocity bins
%       y: mean spike rate / membrane potential
%       w: weights (nr of observations)
%
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  Startpoint for Vm fit = [1 -0.005 0 -50]
%  Startpoint for rate fit = [5 0.05 0 0]


%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell(1,1);
gof = struct( 'sse', cell( 1, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: sigmoidal fit.
[xData, yData, weights] = prepareCurveData(x,y,w);

% Set up fittype and options.
ft = fittype( 'a+b*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
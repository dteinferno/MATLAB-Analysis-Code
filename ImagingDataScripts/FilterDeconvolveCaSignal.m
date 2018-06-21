function dffMatNew = FilterDeconvolveCaSignal(data, dt, tau, filtParam1, filtParam2)
% Filters (sgolay) and deconvolves Ca traces in data.dffMat using Josh
% Vogelstein's deconvolution code (fast_oopsi). See:
% Vogelstein et al., J. Neurophysiol. 2010
%
% Input:
%   data.dffMat: matrix of time series (n * t)
%   dt: time step (sec/frame)
%   filtParam1: sgolayfilt first parameter
%   filtParam2: sgolayfilt second parameter
%
% Output:
%   dffMatNew: filtered, convolved matrix of Ca data
% 

% dt = 0.1175; % seconds/frame
% tau = 0.4/log(2); % t-half/ln(2)
%filtParam1 = 3;
%filtParam2 = 7;

dffMat = sgolayfilt(data.dffMat, filtParam1, filtParam2, [], 2);
nROIs = size(dffMat,1);
dffMatNew = dffMat;

if 0 % Uncomment if you want to use Vogelstein deconvolution algorithm
for i = 1:nROIs
    V.dt = dt;
    P.fast_iter_max = 100;
    P.gam = 1 - dt/tau;
    
    %dfFilt = sgolayfilt(dffMat(i,:), filtParam1, filtParam2);
    [n_best, P_best, V, C] = fast_oopsi(dffMat(i,:), V, P);
    dffMatNew(i,:) = n_best;
end
end

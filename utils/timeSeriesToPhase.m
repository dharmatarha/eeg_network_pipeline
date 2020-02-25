function [phaseData] = timeSeriesToPhase(timeSeriesData)
%% Calculate instantaneous phase of a real-valued signal
% 
% USAGE: phaseData = timeSeriesToPhase(timeSeriesData)
%
% Calculates instantaneous phase of the real valued input singal. 
% It works either with a row vector input (single variable / time series) 
% or a matrix where each row is treated as a variable / time series.
%
% Works by performing Hlbert-transform, and calculating the angle of the
% analytic signal.
%
% Input: 
% timeSeriesData    - Numerical matrix (real), input data, each row is a
%               time series. Phase calculation is done for each row separately.
%
% Output:
% phaseData         - numerical matrix (real) with the same size as input 
%               matrix "timeSeriesData", contains instantaneous phase data.


%% Input checks

if ~ismatrix(timeSeriesData) || ~isreal(timeSeriesData)
    error('Function timeSeriesToPhase expects a numerical matrix of reals as input!');
end

% rows are variables - give a warning if there seems to be more
% variables than time points
if size(timeSeriesData, 1) > size(timeSeriesData, 2)
    warning(['There seems to be more time series than data points in ',...
        'each series - are you sure rows are separate time series?'])
end


%% Instantaneous phase calculation

analyticSignal = hilbert(timeSeriesData');
phaseData = angle(analyticSignal);
phaseData = phaseData';

return
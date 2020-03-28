function [surrConnData] = getSurrConn(realData, surrNo)
%% Generate surrogate connectivity data from real time series data
% 
% USAGE: surrConnData = getSurrConn(realData, surrNo=10000)
%
% Generates surrogate connectivity data from the input vector/matrix using 
% phase scrambling (phase randomization) and calculates connectivity using
% the Phase Locking Value (PLV).
% It works either with row vector input (single variable / time series) or 
% a matrix where each row is treated as a variable / time series.
%
% Input: 
% realData         - Numerical matrix (real), input data, each row is a
%                   time series.
% surrNo           - Number of surrogates. Defaults to 10000.
%
% Output:
% surrConnData     - Numerical matrix (real) with connectivity values
%                   with all possible combinations between time series in 
%                   matrix "realData".
%

%% Input checks

if nargin == 1
    surrNo = 10000;
elseif nargin ~= 2
    error('Function getSurrConn expects input args "realData" and "numberOfSurrogates"!');
end
if ~ismatrix(realData) || ~isreal(realData)
    error('Function getSurrConn expects a numerical matrix of reals as input!');
end
if ~ismember(surrNo, 1:10^5)
    error('Input arg "surrNo" is expected to be between 1:10e5!');
end

% rows are variables - give a warning if there seems to be more
% variables than time points
if size(realData, 1) > size(realData, 2)
    warning(['There seems to be more time series than data points in ',...
        'each series - are you sure rows are separate time series?'])
end


%% Surrogate connectivity calculation

% number of ROIs
roiNo = size(realData, 1);
% variable for storing surrogate connectivity results
surrConnData = nan(surrNo, roiNo, roiNo);

% loop through surrogate data sets: generate, turn to phases, calculate
% connectivity
for surrogateNumber = 1 : surrNo
    yPhaseRand = phaseScramble(realData);
    surrogatePhaseData = timeSeriesToPhase(yPhaseRand);
    surrConnData(surrogateNumber, :, :) = plv(surrogatePhaseData, 0);
end


return
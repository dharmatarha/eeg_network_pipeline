function [surrogateConnectivityData] = generateSurrogateConnectivityData(realData, numberOfSurrogates)
%% Generate surrogate connectivity data from real time series data
% 
% USAGE: surrogateConnectivityData = generateSurrogateConnectivityData(realData, numberOfSurrogates=10000)
%
% Generate surrogate connectivity data of the input matrix using phase scrambing.
% It works either with row vector input (single variable / time series) or a matrix
% where each row is treated as a variable / time series.
%
% Works by performing phase scrambling (phase randomization), and calculates 
% connectivity using the Phase Locking Value (PLV).
%
% Input: 
% realData                      - Numerical matrix (real), input data, each row is a
%                           time series.
% numberOfSurrogates            - Number of surrogates. Defaults to 10000.
%
% Output:
% surrogateConnectivityData     - Numerical matrix (real) with connectivity values
%                          with all possible combinations between time series in matrix 
%                          "realData".
%


%% Input checks

if nargin == 1
    numberOfSurrogates = 10000;
elseif nargin ~= 2
    error('Function generateSurrogateConnectivityData expects input args "realData" and "numberOfSurrogates"!');
end
if ~ismatrix(realData) || ~isreal(realData)
    error('Function generateSurrogateConnectivityData expects a numerical matrix of reals as input!');
end
if ~ismember(numberOfSurrogates, 1:10^5)
    error('Input arg "numberOfSurrogates" is expected to be between 1:10e5!');
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
surrogateConnectivityData = nan(numberOfSurrogates, roiNo, roiNo);

% loop through surrogate data sets: generate, turn to phases, calculate
% connectivity
for surrogateNumber = 1 : numberOfSurrogates
    yPhaseRand = phaseScramble(realData);
    surrogatePhaseData = timeSeriesToPhase(yPhaseRand);
    surrogateConnectivityData(surrogateNumber, :, :) = plv(surrogatePhaseData, 0);
end


return
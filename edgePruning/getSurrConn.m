function [surrConnData] = getSurrConn(realData, varargin)
%% Generate surrogate connectivity data from real time series data
% 
% USAGE: surrConnData = getSurrConn(realData, surrNo=10^4, method='iplv')
%
% Generates surrogate connectivity data from the input matrix using 
% phase scrambling (phase randomization) and calculates connectivity using
% the specified method (one of {'pli', 'plv', 'iplv'}).
% It works either with row vector input (single variable / time series) or 
% a matrix where each row is treated as a variable / time series.
%
% Mandatory input: 
% realData         - Numeric matrix (real) sized (no. of channels/rois X samples). 
%                   Multichannel input data, with each row a separate 
%                   channel / time series.
%
% Optional inputs:
% surrNo           - Numeric value, one of 1:10^5. Defines number of 
%                   surrogates. Defaults to 10^4.
% method           - String, one of {'plv', 'iplv', 'pli'}. Specifies a 
%                   connectivity measure compatible with phase data. 
%                   Defaults to 'iplv'.
%
% Output:
% surrConnData     - 3D numeric array (real) sized 
%                   (no. of surrogates X no. of channels/rois X no. of channels/rois). 
%                   Contains the connectivity values for all possible 
%                   channel pairings for each surrogate data set.
%

%% Input checks

if ~ismember(nargin, 1:3)
    error('Function getSurrConn requires input arg "realData" while args "surrNo" and "method" are optional!');
end
if ~ismatrix(realData) || ~isreal(realData)
    error('Function getSurrConn expects a numerical matrix of reals as input arg "realData"!');
end
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'pli', 'iplv', 'plv'})
            method = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 1:10^5)
            surrNo = varargin{v};
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to args "surrNo" and "method"!']);
        end
    end
else
    surrNo = 10^4;
    method = 'iplv';
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

% loop through surrogate data sets: generate surrogate, turn to phases, calculate
% connectivity
for surrIdx = 1 : surrNo
    yPhaseRand = phaseScramble(realData);  % returns phase-scrambled version of input data
    surrogatePhaseData = timeSeriesToPhase(yPhaseRand);  % extracts instantaneous phase from analyticial signal
    % connectivity measure is specified by input arg "method"
    switch method
        case 'pli'
            surrConnData(surrIdx, :, :) = pli(surrogatePhaseData, 0);
        case 'plv'
            surrConnData(surrIdx, :, :) = plv(surrogatePhaseData, 0);
        case 'iplv'
            surrConnData(surrIdx, :, :) = iplv(surrogatePhaseData, 0);
    end
end


return
function [surrConnData] = getSurrConn(realData, varargin)
%% Generate surrogate connectivity data from real time series data
% 
% USAGE: surrConnData = getSurrConn(realData, surrNo=10^4, method='iplv', lpFilter=[])
%
% Generates surrogate connectivity data from the input matrix using 
% phase scrambling (phase randomization) and calculates connectivity using
% the specified method 
% (one of {'pli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}).
%
% It works either with row vector input (single variable / time series) or 
% a matrix where each row is treated as a variable / time series.
%
%
% Mandatory input: 
% realData         - Numeric matrix (real) sized (no. of channels/rois X samples). 
%                   Multichannel input data, with each row a separate 
%                   channel / time series.
%
% Optional inputs:
% surrNo            - Numeric value, one of 1:10^5. Defines number of 
%                   surrogates. Defaults to 10^4.
% method            - String, one of {'pli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}. 
%                   Specifies the connectivity measure. Defaults to 'iplv'.
% lpFilter          -Digital filter object as returned by e.g. designfilt
%                   (part of Signal Processing Toolbox!). Used if "method" 
%                   is either 'ampCorr' or 'orthAmpCorr' 
%                   (envelope correlations) and is for lowpass filtering 
%                   the amplitude envelopes before calculating correlations.
%                   Applied via the Signal Processing Toolbox (!) version of 
%                   the filter command. Defaults to [], meaning no filtering.
%
% Output:
% surrConnData      - 3D numeric array (real) sized 
%                   (no. of surrogates X no. of channels/rois X no. of channels/rois). 
%                   Contains the connectivity values for all possible 
%                   channel pairings for each surrogate data set.
%

%% Input checks

if ~ismember(nargin, 1:4)
    error('Function getSurrConn requires input arg "realData" while args "surrNo", "method" and "lpFilter" are optional!');
end
if ~ismatrix(realData) || ~isreal(realData)
    error('Function getSurrConn expects a numerical matrix of reals as input arg "realData"!');
end
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismember(varargin{v}, 1:10^5) && ~exist('surrNo', 'var')
            surrNo = varargin{v};        
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'pli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}) && ~exist('method', 'var')
            method = varargin{v};
        elseif isa(varargin{v}, 'digitalFilter') && ~exist('lpFilter', 'var')
            lpFilter = varargin{v};
        else
            error(['At least one input arg could not be mapped ',...
                'nicely to args "surrNo", "method" or "lpFilter"!']);
        end
    end
end
if ~exist('surrNo', 'var')
    surrNo = 10^4;
end
if ~exist('method', 'var')
    method = 'iplv';
end
if ~exist('lpFilter', 'var')
    lpFilter = [];
end

% rows are variables - give a warning if there seems to be more
% variables than time points
if size(realData, 1) > size(realData, 2)
    warning(['There seems to be more time series than data points in ',...
        'each series - are you sure rows are separate time series?'])
end

% give a warning if a lowpass filter is supplied with a phase-based
% connectivity measure
if ismember(method, {'pli', 'plv', 'iplv'}) && ~isempty(lpFilter)
    warning(['A lowpass filter object was supplied as input arg ',...
        '(mapped to "lpFilter") but the selected connectivity method is phase-based.',... 
        char(10), 'No filtering will take place.']);
end


%% Surrogate connectivity calculation

% number of ROIs
roiNo = size(realData, 1);
% variable for storing surrogate connectivity results
surrConnData = nan(surrNo, roiNo, roiNo);

% loop through surrogate data sets: generate surrogate, turn to phases if 
% necessary, calculate connectivity
for surrIdx = 1 : surrNo
    surrogateData = phaseScramble(realData);  % returns phase-scrambled version of input data
    
    % for phase-based methods, extract instantaneous phase from analytical signal
    if ismember(method, {'pli', 'plv', 'iplv'})
        surrogatePhaseData = timeSeriesToPhase(surrogateData);
    end
    
    % connectivity measure is specified by input arg "method"
    switch method
        case 'pli'
            surrConnData(surrIdx, :, :) = pli(surrogatePhaseData, 0);
        case 'plv'
            surrConnData(surrIdx, :, :) = plv(surrogatePhaseData, 0);
        case 'iplv'
            surrConnData(surrIdx, :, :) = iplv(surrogatePhaseData, 0);
        case 'ampCorr'
            if ~isempty(lpFilter)
                surrConnData(surrIdx, :, :) = ampCorr(surrogateData, lpFilter, 0);      
            else
                surrConnData(surrIdx, :, :) = ampCorr(surrogateData, 0);
            end
        case 'orthAmpCorr'
            if ~isempty(lpFilter)
                surrConnData(surrIdx, :, :) = orthAmpCorr(surrogateData, lpFilter, 0);
            else
                surrConnData(surrIdx, :, :) = orthAmpCorr(surrogateData, 0);
            end
    end
end


return
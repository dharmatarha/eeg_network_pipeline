function [surrConnData] = getSurrConn(realData, varargin)
%% Generate surrogate connectivity data from real time series data
% 
% USAGE: surrConnData = getSurrConn(realData, surrNo=10^4, method='iplv', lpFilter=[])
%
% Generates surrogate connectivity data from the input matrix using 
% phase scrambling (phase randomization) and calculates connectivity using
% the specified method or methods
% (one or more of {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'}).
%
% It works either with row vector input (single variable / time series) or 
% a matrix where each row is treated as a variable / time series.
%
% Optional inputs are inferred from types and values.
%
% Mandatory input: 
% realData         - Numeric matrix (real) sized (no. of channels/rois X samples). 
%                   Multichannel input data, with each row a separate 
%                   channel / time series.
%
% Optional inputs:
% surrNo            - Numeric value, one of 1:10^5. Defines number of 
%                   surrogates. Defaults to 10^4.
% method            - Either a char array, one of 
%                   {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'}, or 
%                   a cell array of char array. Specifies one or more 
%                   connectivity measures to calculate on the surrogate 
%                   data. Defaults to 'iplv'.
% lpFilter          -Digital filter object as returned by e.g. designfilt
%                   (part of Signal Processing Toolbox!). Used if "method" 
%                   is either 'ampCorr' or 'orthAmpCorr' 
%                   (envelope correlations) and is for lowpass filtering 
%                   the amplitude envelopes before calculating correlations.
%                   Applied via the Signal Processing Toolbox (!) version of 
%                   the filter command. Defaults to [], meaning no filtering.
%
% Output:
% surrConnData      - 4D numeric array (real) sized 
%                   [no. of methods X no. of surrogates X no. of channels/rois X no. of channels/rois]. 
%                   Contains the connectivity values for all connectivity 
%                   methods and all possible channel/roi pairings for each 
%                   surrogate data set.
%


%% Input checks

% check no. of inputs
if ~ismember(nargin, 1:4)
    error('Function getSurrConn requires input arg "realData" while args "surrNo", "method" and "lpFilter" are optional!');
end
% check mandatory input
if ~ismatrix(realData) || ~isreal(realData)
    error('Input arg "realData" should be a numeric matrix of reals!');
end
% check optional inputs
% remember, "method" could be a char array or a cell array of char arrays
% as well
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismember(varargin{v}, 1:10^5) && ~exist('surrNo', 'var')
            surrNo = varargin{v};        
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'}) && ~exist('method', 'var')
            method = varargin{v};
        elseif iscell(varargin{v}) && all(ismember(varargin{v}, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'})) && ~exist('method', 'var')
            method = varargin{v};    
        elseif isa(varargin{v}, 'digitalFilter') && ~exist('lpFilter', 'var')
            lpFilter = varargin{v};
        else
            error(['At least one input arg could not be mapped ',...
                'nicely to args "surrNo", "method" or "lpFilter"!']);
        end
    end
end
% assign default values where necessary
if ~exist('surrNo', 'var')
    surrNo = 10^4;
end
if ~exist('method', 'var')
    method = 'iplv';
end
if ~exist('lpFilter', 'var')
    lpFilter = [];
end

% make sure that "method" is a cell - transform if char array, makes it
% easier to treat it if the type is consistent
if ischar(method)
    method = {method};
end

% rows are variables - give a warning if there seems to be more
% variables than time points
if size(realData, 1) > size(realData, 2)
    warning(['There seems to be more time series than data points in ',...
        'each series - are you sure rows are separate time series?'])
end

% give a warning if a lowpass filter is supplied with only phase-based
% connectivity methods
if all(ismember(method, {'plv', 'iplv', 'ciplv'})) && ~isempty(lpFilter)
    warning(['A lowpass filter object was supplied as input arg ',...
        '(mapped to "lpFilter") but the selected connectivity method(s) is (are) phase-based.',... 
        char(10), 'No filtering will take place.']);
end


%% Surrogate connectivity calculation

% number  of connectivity methods
methodNo = length(method);
% number of ROIs
roiNo = size(realData, 1);
% variable for storing surrogate connectivity results
surrConnData = nan(methodNo, surrNo, roiNo, roiNo);

% loop through surrogate data sets: generate surrogate, turn to phases if 
% necessary, calculate connectivity 
for surrIdx = 1:surrNo
    
    % get phase-scrambled version of input data
    surrogateData = phaseScramble(realData);
    
    % loop through connectivity methods requested
    for methodIdx = 1:methodNo
        % char array of current method
        currentMethod = method{methodIdx};

        % connectivity measure is specified by input arg "method"
        switch currentMethod
            case 'plv'
                surrConnData(methodIdx, surrIdx, :, :) = plv(surrogatePhaseData);
            case 'iplv'
                surrConnData(methodIdx, surrIdx, :, :) = iplv(surrogatePhaseData);
            case 'ciplv'
                surrConnData(methodIdx, surrIdx, :, :) = ciplv(surrogatePhaseData);                
            case 'ampCorr'
                if ~isempty(lpFilter)
                    surrConnData(methodIdx, surrIdx, :, :) = ampCorr(surrogateData, lpFilter);      
                else
                    surrConnData(methodIdx, surrIdx, :, :) = ampCorr(surrogateData);
                end
            case 'orthAmpCorr'
                if ~isempty(lpFilter)
                    surrConnData(methodIdx, surrIdx, :, :) = orthAmpCorr(surrogateData, lpFilter);
                else
                    surrConnData(methodIdx, surrIdx, :, :) = orthAmpCorr(surrogateData);
                end
        end  % switch method
        
    end  % for methodIdx

end  % for surrIdx


return
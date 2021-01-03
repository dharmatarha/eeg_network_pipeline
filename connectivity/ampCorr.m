function res = ampCorr(epochData, varargin)
%% Amplitude / Envelope correlation connectivity measure
%
% USAGE: res = ampCorr(epochData, lpFilter=[], verbosity=1)
%
% Calculates envelope correlations on epoch-level data. Expects 
% real-valued data as input in ROIs/channels X samples format. 
%
% As correlation is symmetric, only the upper triangle of the output matrix
% is populated.
%
% Mandatory input:
% epochData     - Numeric matrix where each row is a separate channel / ROI 
%               time series and columns correspond to samples. 
%
% Optional inputs:
% lpFilter      - Digital filter object as returned by e.g. designfilt
%               (part of Signal Processing Toolbox!). Used for lowpass
%               filtering the amplitude envelopes before correlations.
%               Applied via the Signal Processing Toolbox (!) version of 
%               the filter command. Defaults to [], meaning no filtering. 
% verbosity     - Verbosity. If 1, it prints to command window, 0 means
%               silence. Default is 1.
%
% Output:
% res           - Numeric matrix of connectivity values, sized
%               ROIs/channels X ROIs/channels, so that entry i,j is the
%               envelope correlation of ROIs/channels i and j. Only the
%               upper triangle above the main diagonal is populated with 
%               values, rest is NaN.
%


%% Input checks

if ~ismember(nargin, 1:3)
    error('Function ampCorr requires input arg "epochData" while args "lpFilter" and "v" are optional!');
end
if ~isnumeric(epochData) || ~ismatrix(epochData)
    error('Input arg "epochData" should be  numeric matrix!');
end
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isa(varargin{v}, 'digitalFilter') && ~exist('lpFilter', 'var')
            lpFilter = varargin{v};
        elseif ismember(varargin{v} ,[0 1]) && ~exist('verbosity', 'var')
            verbosity = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "lpFilter" or "verbosity"!');            
        end
    end
end
if ~exist('lpFilter', 'var')
    lpFilter = [];
end
if ~exist('verbosity', 'var')
    verbosity = 1;
end

% check input data size
[roiNo, sampleNo] = size(epochData);
if roiNo >= sampleNo
    warning('Input data seems to have many ROIs/channels relative to the number of samples. We proceed but it is suspicious.');
end

% user message
if verbosity
    disp([char(10), 'Called ampCorr on data with ', num2str(roiNo),... 
        ' ROIs/channels, each with ', num2str(sampleNo), ' samples']);
end


%% Get envelope, filter, get correlations

% envelope based on analytical signal
epochDataEnv = envelope(epochData(:, :)')';  % we keep the dim order channels/ROIs X samples

% filter envelopes if there was a filter object supplied
if ~isempty(lpFilter)
    epochDataEnv = filter(lpFilter, epochDataEnv')';  % we keep the dim order channels/ROIs X samples
end

% correlations
res = corr(epochDataEnv')';  % we keep the dim order channels/ROIs X samples

% lower triangle is set to NaN, as with symmetric phase-based measures  
res(tril(true(roiNo))) = NaN; 

% user message
if verbosity
    disp(['Calculated ampCorr for ', num2str(roiNo*(roiNo-1)/2), ' ROI/channel pairings']);
end


return
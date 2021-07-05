function res = ampCorr(data, varargin)
%% Amplitude / Envelope correlation connectivity measure
%
% USAGE: res = ampCorr(data, lpFilter=[])
%
% Calculates envelope correlations on epoch-level data. Expects 
% real-valued data as input with dimensions ROIs/channels X samples. 
%
% As correlation is symmetric, only the upper triangle of the output matrix
% is populated.
%
% Mandatory input(s):
% data          - Numeric matrix, real valued. Its dimensions are channels
%               / ROI time series X samples. 
%
% Optional input(s):
% lpFilter      - Digital filter object as returned by e.g. designfilt
%               (part of Signal Processing Toolbox!). Used for lowpass
%               filtering the amplitude envelopes before correlations.
%               Applied via the Signal Processing Toolbox (!) version of 
%               the filter command. Defaults to [], meaning no filtering. 
%
% Output:
% res           - Numeric matrix of connectivity values, sized
%               channels X channels, so that entry i,j is the
%               envelope correlation of ROIs/channels i and j. Only the
%               upper triangle above the main diagonal is populated with 
%               values, rest is NaN.
%


%% Input checks

% number of input args
if ~ismember(nargin, 1:2)
    error('Function ampCorr requires input arg "data" while arg "lpFilter" is optional!');
end
% check mandatory input
if ~isnumeric(data) || ~isreal(data) || ~ismatrix(data)
    error('Input arg "data" should be numeric matrix of reals!');
end
% check / assign optional input
if nargin == 2
    if isa(varargin{1}, 'digitalFilter') || isempty(varargin{1})
        lpFilter = varargin{1};
    else
        error('Optional arg "lpFilter" should be either a digital filter object or empty!');            
    end
else
    lpFilter = [];
end

% sanity check on input data size
[channelNo, sampleNo] = size(data);
if channelNo >= sampleNo
    warning('Input data has equal or larger number of channels than samples. We proceed but it is suspicious!');
end


%% Get envelope, filter, get correlations

% envelope based on analytical signal
dataEnv = envelope(data(:, :)')';  % we keep the dim order channels/ROIs X samples

% filter envelopes if there was a filter object supplied
if ~isempty(lpFilter)
    dataEnv = filter(lpFilter, dataEnv')';  % we keep the dim order channels X samples
end

% correlations
res = corr(dataEnv')';  % we keep the dim order channels X samples

% lower triangle is set to NaN, as with symmetric phase-based measures  
res(tril(true(channelNo))) = NaN; 


return
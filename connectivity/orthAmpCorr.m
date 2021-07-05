function res = orthAmpCorr(data, varargin)
%% Pairwise-orthogonalized envelope correlation
%
% USAGE: res = orthAmpCorr(data, lpFilter=[])
%
% Calculates envelope correlations on pairwise-orthogonalized ROI/channel
% data. Expects real-valued data as input with dimensions ROIs/channels 
% X samples. 
%
% As pairwise-orthogonalization is noncommutative, it is performed in both
% directions and the average of the two connectivity values is returned in
% the output arg.
%
% Mandatory input(s):
% data          - Numeric matrix, real valued. Its dimensions are channels
%               / ROI time series X samples.
%
% Optional inputs:
% lpFilter      - Digital filter object as returned by e.g. designfilt
%               (part of Signal Processing Toolbox!). Used for lowpass
%               filtering the amplitude envelopes before correlations.
%               Applied via the Signal Processing Toolbox (!) version of 
%               the filter command. Defaults to [], meaning no filtering. 
%
% Output:
% res           - Numeric matrix of connectivity values, sized
%               channels X channels, so that entry i,j is the
%               envelope correlation of channels i and j. Only the
%               upper triangle above the main diagonal is populated with 
%               values, rest is NaN.
%
%
% Relevant papers:
% Coquelet et al., 2020. Comparing MEG and high-density EEG for intrinsic
% functional connectivity mapping. NeuroImage.
%


%% Input checks

% number of input args
if ~ismember(nargin, 1:2)
    error('Function orthAmpCorr requires input arg "data" while arg "lpFilter" is optional!');
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


%% Loops through all ROI/channel pairings

% transpose to samples X channels
data = data';

% temp var for connectivity res
res = nan(channelNo);

% vector norms (Euclidean) of all channels 
channelNorms = sqrt(sum(data.^2, 1));  % channelNorms is a row vector

% get envelopes for all channels (across samples)
dataEnv = envelope(data);  % dimenions here are samples X channels

% filter all channels if lowpass filter was provided
if ~isempty(lpFilter)
    dataEnv = (filter(lpFilter, dataEnv));
end

% loops trough ROIs
for channelIdx = 1:channelNo 
    
    % select one channel to correlate with all others, replicate it
    % channelNo times
    channelData = repmat(data(:, channelIdx), [1, 62]);  

    % scalar projections of current channel data on all channels
    projScalars = dot(channelData, data)./channelNorms;  % projScalar is 1 X channelNo
    
    % projection vectors of current channel data on all channels
    projVectors = (data./channelNorms).*projScalars;  
    
    % orthogonalized current channel with regards to all other channels
    orthVectors = channelData - projVectors;
    
    % get envelopes
    channelDataEnv = envelope(orthVectors);  
    
    % filter envelopes
    if ~isempty(lpFilter)
        channelDataEnv = (filter(lpFilter, channelDataEnv));
    end
    
    % get correlations
    rhos = corr(channelDataEnv, dataEnv);
    
    % take only the diagonal for the results matrix
    res(channelIdx, :) = diag(rhos);
    
end

% average the upper and lower triangles into the upper
% triangle
res = (triu(res, 1) + tril(res, -1)')/2;
% lower triangle is set to NaN
res(tril(true(channelNo))) = NaN;  


return








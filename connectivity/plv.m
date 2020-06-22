function plvRes = plv(epochData, v)
%% Phase-Locking value
%
% USAGE: plvRes = plv(epochData, v=1)
%
% Function to calculate phase-locking values (PLVs) between a set of 
% channels / time series. IMPORTANT: Works on phase values!
%
% Input(s):
% epochData     - The input is a matrix of phase values (numeric between 
%               -1 +1 pi) where each row is a separate channel / time 
%               series and columns correspond to samples.  
% v             - Verbosity. If 1, it prints to command window, 0 means
%               silence. Default is 1.
%
% Output(s):
% plvRes        - Matrix (sized no. of channels X no. of channels) of PLV
%               values where entry i,j is PLV between channels i and j. As
%               the matrix is symmetric across the diagonal, there are
%               values only in the upper triangle, the rest is NaN.
%
% Relevant papers:
% Lachaux et al., 1999. Measuring phase synchrony in brain signals. 
%   Hum. Brain Mapp.
% Mormann et al., 2000. Mean phase coherence as a measure for phase 
%   synchronization and its application to the EEG of epilepsy patients. 
%   Physica D.
%


%% Input checks

% number of input args
if nargin == 1 
    v = 1;
elseif nargin ~= 2
    error('Function plv requires input arg "epochData" and optional arg "v"!');
end
% check verbosity
if ~ismembertol(v, [0 1])
    error('Input arg "v" is either 0 or 1!');
end
% check input data dimensionality
if length(size(epochData)) ~= 2
    error('Input arg "epochData" is expected to be a 2D matrix (channels/ROIs X samples)!');
end
% are input values between -pi +pi ?
if any(any(epochData > pi)) || any(any(epochData < -pi))
    error('Elements of input arg "epochData" are >pi or <-pi. Are you sure these are phase values?');
end

% get number of channels and samples
channelNo = size(epochData, 1);
sampleNo = size(epochData, 2);

% user message
if v
    disp([char(10), 'Called plv on data with ', num2str(channelNo),... 
        ' channels, each with ', num2str(sampleNo), ' samples']);
end


%% Loop across channels

% preallocate results variable
plvRes = nan(channelNo, channelNo);

% go through channels/ROIs
for channelOne = 1:channelNo-1
    % fill a matrix with repetitions of the data from current channel/ROI
    d1 = repmat(epochData(channelOne, :), [channelNo-channelOne, 1]);
    % calculate plv on matrices (= multiple channel pairings)
    plvRes(channelOne, channelOne+1:end) = abs(sum(exp(1i*(d1-epochData(channelOne+1:end, :))), 2)/sampleNo);
    
end

% user message
if v
    disp(['Calculated plv for ', num2str(channelNo*(channelNo-1)/2), ' channel pairings']);
end

return


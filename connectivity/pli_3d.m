function pliRes = pli_3d(epochData, v)
%% Phase-Locking Index for 3D numeric array as input
%
% USAGE: pliRes = pli_3d(epochData, v=1)
%
% Function to calculate phase-locking indices (PLIs) between a set of 
% channels / time series across epochs / trials. That is, the current
% function can estimate connectivity separately for samples (time points).
%
% IMPORTANT: Works on phase values!
%
% Input(s):
% epochData     - The input is a 3D array of phase values (numeric between 
%               -1 +1 pi) with dimensions 
%               [epochs/trials X channels/ROIs X samples]. PLI is 
%               estimated across epochs/trials, separetely for samples.  
% v             - Verbosity. If 1, it prints to command window, 0 means
%               silence. Default is 1.
%
% Output(s):
% pliRes        - 3D numeric array with dimensions 
%               [channels/ROIs X channels/ROIs X samples]. Contains PLI
%               values where entry i, j, k is PLI between channels i and j 
%               for sample k. For each sample ("k"), the connectivity 
%               matrix is symmetric across the diagonal and there are
%               values only in the upper triangle, the rest is NaN.
%
% Relevant papers:
% Stam et al., 2007. Phase lag index: Assessment of functional connectivity 
%   from multi channel EEG and MEG with diminished bias from common sources. 
%   Hum. Brain Mapp.
% See also:
% Nolte et al., 2004. Identifying true brain interaction from EEG data 
%   using the imaginary part of coherency. Clin. Neurophys.
%


%% Input checks

% number of input args
if nargin == 1 
    v = 1;
elseif nargin ~= 2
    error('Function pli_3d requires input arg "epochData" and optional arg "v"!');
end
% check mandatory input
if ~isnumeric(epochData) || numel(size(epochData)) ~= 3
   error('Input arg epochData should be a 3D numeric array!'); 
end
% check optional arg
if ~ismembertol(v, [0 1])
    error('Input arg "v" is either 0 or 1!');
end
% are input values between -pi +pi ?
if any(epochData(:) > pi) || any(epochData(:) < -pi)
    error('Elements of input arg "epochData" are >pi or <-pi. Are you sure these are phase values?');
end

% get number of epochs, channels and samples
[epochNo, channelNo, sampleNo] = size(epochData);

% user message
if v
    disp([char(10), 'Called pli_3d on data with ', num2str(channelNo),... 
        ' channels, each with ', num2str(sampleNo), ' samples']);
end


%% Loop across channel pairings

% preallocate results variable
pliRes = nan(channelNo, channelNo, sampleNo);


%% Calculations, for-loop version 
%
% % Left here as a ground-truth implementation to compare any faster
% % solution with
%
% % loop through samples
% for sampleIdx = 1:sampleNo
%     % go through channels/ROIs
%     for channelOne = 1:channelNo-1
%         % fill a matrix with repetitions of the data from current channel/ROI
%         d1 = repmat(squeeze(epochData(:, channelOne, sampleIdx))', [channelNo-channelOne, 1]);
%         % calculate pli on matrices (= multiple channel pairings)
%         pliRes(channelOne, channelOne+1:end, sampleIdx) = abs(mean(sign(d1-squeeze(epochData(:, channelOne+1:end, sampleIdx))'), 2));
%     end  % for channelOne
% end  % for sampleIdx


%% Calculations, avoiding loop over samples

% go through channels/ROIs
for channelOne = 1:channelNo-1
    % fill a 3D array with repetitions of the data from current channel/ROI
    d1 = repmat(epochData(:, channelOne, :), [1, channelNo-channelOne, 1]);
    % calculate pli on whole arrays (= multiple channel pairings)
    pliRes(channelOne, channelOne+1:end, :) = abs(mean(sign(d1 - epochData(:, channelOne+1:end, :)), 1));
end  % for channelOne

% user message
if v
    disp(['Calculated pli for ', num2str(channelNo*(channelNo-1)/2), ' channel pairings']);
end


return


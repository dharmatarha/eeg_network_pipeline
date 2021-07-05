function plvRes = plv(data)
%% Phase-Locking value
%
% USAGE: plvRes = plv(data)
%
% Function to calculate phase-locking values (PLVs) across a set of 
% channels / ROI time series. Supports PLV calculation in time series, 
% not across trials.
%
% Returns PLV values only in the upper triangle.
%
% IMPORTANT: Expects real-valued data, unlike earlier versions!
% 
% Mandatory input(s):
% data          - Numeric matrix, real valued. Its dimensions are channels
%               X samples.
%
% Output(s):
% plvRes        - Numeric matrix, real valued. Contains PLV values. Its 
%               dimensions are channels X channels). Entry i,j is PLV 
%               between channels i and j.
%
% Relevant papers:
% Lachaux et al., 1999. Measuring phase synchrony in brain signals. 
%   Hum. Brain Mapp.
% Mormann et al., 2000. Mean phase coherence as a measure for phase 
%   synchronization and its application to the EEG of epilepsy patients. 
%   Physica D.
% Bruna et al., 2018. Phase locking value revisited: teaching new tricks 
%   to an old dog. J. Neural Eng. 
%
% NOTES:
% Previous implementation was based on Mormann et al. (2000), currently
% implements the much faster method derived in Bruna et al. (2018). 
%


%% Input checks

% number of input args
if nargin ~= 1
    error('Function plv requires input arg "data"!');
end
% check mandatory input
if ~isnumeric(data) || ~isreal(data) || ~ismatrix(data)
    error('Input arg "data" should be a 2D numeric matrix of reals (channels X samples)!');
end

% sanity check on input data size
[channelNo, sampleNo] = size(data);
if channelNo >= sampleNo
    warning('Input data has equal or larger number of channels than samples. We proceed but it is suspicious!');
end


%% Calculation

% NEW, FAST METHOD BASED ON BRUNA ET AL. (2018):

% get analytic signal, across samples
dataAnalytic = hilbert(data')';  % transposes keep the dims as channels X samples
% normalize with magnitude
normedData = dataAnalytic ./ abs(dataAnalytic);
% plv as matrix multiplication
plvRes = abs(normedData * normedData') / sampleNo;

% lower triangle is set to NaN, as with symmetric phase-based measures  
plvRes(tril(true(channelNo))) = NaN; 


% % OLD, SLOWER METHOD BASED ON MORMANN ET AL. (2000):
% 
% % get analytic signal, across samples
% dataAnalytic = hilbert(data')';  % transposes keep the dims as channels X samples
% dataPhase = angle(dataAnalytic);  % phase data
% 
% % preallocate results variable
% plvRes = nan(channelNo, channelNo);
% 
% % go through channels/ROIs
% for currentChannel = 1:channelNo-1
%     % fill a matrix with repetitions of the data from current channel/ROI
%     d1 = repmat(dataPhase(currentChannel, :), [channelNo-currentChannel, 1]);
%     % calculate plv on matrices (= multiple channel pairings)
%     plvRes(currentChannel, currentChannel+1:end) = abs(sum(exp(1i*(d1-dataPhase(currentChannel+1:end, :))), 2)/sampleNo);
%     
% end



return


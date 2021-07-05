function iplvRes = iplv(data)
%% Imaginary part of the Phase-Locking Value (iPLV)
%
% USAGE: iplvRes = iplv(data)
%
% Function to calculate imaginary part of phase-locking values (iPLVs) 
% across a set of channels / ROI time series. Supports iPLV calculation in 
% time series, not across trials.
%
% Returns iPLV values only in the upper triangle.
%
% IMPORTANT: Expects real-valued data, unlike earlier versions!
% 
% Mandatory input(s):
% data          - Numeric matrix, real valued. Its dimensions are channels
%               X samples.
%
% Output(s):
% iplvRes       - Numeric matrix, real valued. Contains iPLV values. Its 
%               dimensions are channels X channels). Entry i,j is iPLV 
%               between channels i and j.
%
% Relevant papers:
% Palva, S., & Palva, J. M., 2012. Discovering oscillatory interaction 
%   networks with M/EEG: challenges and breakthroughs. Trends in cog. sci.
% Palva et al., 2018. Ghost interactions in MEG/EEG source space: A note 
%   of caution on inter-areal coupling measures. Neuroimage.
% Bruna et al., 2018. Phase locking value revisited: teaching new tricks 
%   to an old dog. J. Neural Eng. 
% Check also the papers mentioned in the help of plv.m for the background
% on regular PLV.
%
% NOTES:
% Previous implementation was based on Palva et al. (2012, 2018),
% currently implements the much faster method derived in Bruna et al.
% (2018). 
%


%% Input checks

% number of input args
if nargin ~= 1
    error('Function iplv requires input arg "data"!');
end
% check mandatory input
if ~isnumeric(data) || ~isreal(data) || length(size(data)) ~= 2
    error('Input arg "data" should be a 2D numeric matrix of reals (channels X samples)!');
end

% sanity check on input data size
[channelNo, sampleNo] = size(data);
if channelNo >= sampleNo
    warning('Input data has equal or larger number of channels than samples. We proceed but it is suspicious!');
end


%% Calculations

% NEW, FAST METHOD BASED ON BRUNA ET AL. (2018):

% get analytic signal, across samples
dataAnalytic = hilbert(data')';  % transposes keep the dims as channels X samples
% normalize with magnitude
normedData = dataAnalytic ./ abs(dataAnalytic);
% iplv as matrix multiplication
iplvRes = abs(imag(normedData * normedData') / sampleNo);

% lower triangle is set to NaN, as with symmetric phase-based measures  
iplvRes(tril(true(channelNo))) = NaN; 


% % OLD, SLOWER METHOD BASED ON MORMANN ET AL. (2000):
% 
% % get analytic signal, across samples
% dataAnalytic = hilbert(data')';  % transposes keep the dims as channels X samples
% dataPhase = angle(dataAnalytic);  % phase data
% % preallocate results variable
% iplvRes = nan(channelNo, channelNo);
% 
% % go through channels/ROIs
% for channelOne = 1:channelNo-1
%     % fill a matrix with repetitions of the data from current channel/ROI
%     d1 = repmat(dataPhase(channelOne, :), [channelNo-channelOne, 1]);
%     % calculate iplv on matrices (= multiple channel pairings)
%     iplvRes(channelOne, channelOne+1:end) = abs(imag(sum(exp(1i*(d1-dataPhase(channelOne+1:end, :))), 2)/sampleNo));
%     
% end


return


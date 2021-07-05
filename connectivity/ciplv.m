function ciplvRes = ciplv(data)
%% Corrected Imaginary part of the Phase-Locking Value (ciPLV)
%
% USAGE: ciplvRes = ciplv(data)
%
% Function to calculate the corrected imaginary part of phase-locking 
% value (ciPLV) across a set of channels / ROI time series. Supports ciPLV 
% calculation in time series, not across trials.
% 
% Returns ciPLV values only in the upper triangle.
%
% IMPORTANT: Expects real-valued data, unlike earlier versions!
% 
% Mandatory input(s):
% data          - Numeric matrix, real valued. Its dimensions are channels
%               X samples.
%
% Output(s):
% ciplvRes      - Numeric matrix, real valued. Contains ciPLV values. Its 
%               dimensions are channels X channels). Entry i,j is ciPLV 
%               between channels i and j.
%
% Relevant papers:
% Bruna et al., 2018. Phase locking value revisited: teaching new tricks 
%   to an old dog. J. Neural Eng. 
% Check also the papers mentioned in the help of plv.m and iplv.m for 
% the background on regular PLV and on iPLV.
% 
% NOTE: 
% - The magnitude of the imaginery part (nominator) and of the real
% part (in denominator) is taken as in the implementation by MNE-TOOLS. 
% Otherwise the scale is -1 : +1, while the measure is undirected....  
% - Intuitively, the ciPLV corrects the coupling magnitude closer to 1 if 
% there is consistent coupling, but the mean phase difference is close to 0. 
% It only penalizes heavily truly 0 / 180 degree phase-difference coupling.


%% Input checks

% number of input args
if nargin ~= 1
    error('Function ciplv requires input arg "data"!');
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


%% Calculation

% get analytic signal, across samples
dataAnalytic = hilbert(data')';  % transposes keep the dims as channels X samples
% normalize with magnitude
normedData = dataAnalytic ./ abs(dataAnalytic);
% channel-pairing-wise dot product of normalized analytic data, with matrix
% multiplication
tmp = (normedData * normedData');
% corrections: extract the magnitude of the imaginary parts and normalize
% with the magnitude of the real part
ciplvRes = abs(imag(tmp)/sampleNo) ./ (sqrt(1-(abs(real(tmp))/sampleNo).^2));

% lower triangle is set to NaN, as with symmetric phase-based measures  
ciplvRes(tril(true(channelNo))) = NaN; 

return

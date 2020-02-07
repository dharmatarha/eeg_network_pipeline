function plvRes = plv(epochData)
%% Phase-Locking value
%
% USAGE: plvRes = plv(epochData)
%
% Function to calculate phase-locking values (PLVs) between a set of 
% channels / time series. IMPORTANT: Works on phase values!
%
% Input(s):
% epochData     - The input is a matrix of phase values (numeric between 
%               -1 +1 pi) where each row is a separate channel / time 
%               series and columns correspond to samples.  
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
if nargin ~= 1
    error('Function plv requires input arg "epochData"!');
end
% are input values between -pi +pi ?
if any(any(epochData > pi)) || any(any(epochData < -pi))
    error('Elements of input arg "epochData" are >pi or <-pi. Are you sure these are phase values?');
end

% get number of channels and samples
channelNo = size(epochData, 1);
sampleNo = size(epochData, 2);

% user message
disp([char(10), 'Called plv on data with ', num2str(channelNo),... 
    ' channels, each with ', num2str(sampleNo), ' samples']);


%% Loop across channel pairings

% preallocate results variable
plvRes = nan(channelNo, channelNo);
% variable for storing channel pairings we already calculated PLV for
pastPairings = nan(channelNo*(channelNo-1)/2, 2);
% channel pairing counter
counter = 0;

% loops over channels
for channelOne = 1:channelNo
    for channelTwo = 1:channelNo
        
        % only calculate PLV if channel numbers do not match and have not
        % been encountered before
        if channelOne ~= channelTwo && ~ismember([channelOne, channelTwo], pastPairings, 'rows') && ~ismember([channelTwo, channelOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [channelOne, channelTwo];
            
            % We use the Mormann et al. version here:
            plvRes(channelOne, channelTwo) = abs(sum(exp(1i*(epochData(channelOne, :)-epochData(channelTwo, :))))/sampleNo);            
            
        end  % pairings if-then
        
    end  % channelTwo for loop  
end  % channelOne for loop

% user message
disp(['Calculated plv for ', num2str(channelNo*(channelNo-1)/2), ' channel pairings']);


return


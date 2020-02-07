function pliRes = pli(epochData)
%% Phase-Locking Index
%
% USAGE: pliRes = pli(epochData)
%
% Function to calculate phase-locking indices (PLIs) between a set of 
% channels / time series. IMPORTANT: Works on phase values!
%
% Input(s):
% epochData     - The input is a matrix of phase values (numeric between 
%               -1 +1 pi) where each row is a separate channel / time 
%               series and columns correspond to samples.  
%
% Output(s):
% pliRes        - Matrix (sized no. of channels X no. of channels) of PLI
%               values where entry i,j is PLI between channels i and j. As
%               the matrix is symmetric across the diagonal, there are
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
if nargin ~= 1
    error('Function pli requires input arg "epochData"!');
end
% are input values between -pi +pi ?
if any(any(epochData > pi)) || any(any(epochData < -pi))
    error('Elements of input arg "epochData" are >pi or <-pi. Are you sure these are phase values?');
end

% get number of channels and samples
channelNo = size(epochData, 1);
sampleNo = size(epochData, 2);

% user message
disp([char(10), 'Called pli on data with ', num2str(channelNo),... 
    ' channels, each with ', num2str(sampleNo), ' samples']);


%% Loop across channel pairings

% preallocate results variable
pliRes = nan(channelNo, channelNo);
% variable for storing channel pairings we already calculated PLI for
pastPairings = nan(channelNo*(channelNo-1)/2, 2);
% channel pairing counter
counter = 0;

% loops over channels
for channelOne = 1:channelNo
    for channelTwo = 1:channelNo
        
        % only calculate PLI if channel numbers do not match and have not
        % been encountered before
        if channelOne ~= channelTwo && ~ismember([channelOne, channelTwo], pastPairings, 'rows') && ~ismember([channelTwo, channelOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [channelOne, channelTwo];
            
            % We use the Mormann et al. version here:
            pliRes(channelOne, channelTwo) = abs(mean(sign(epochData(channelOne, :)-epochData(channelTwo, :))));            
            
        end  % pairings if-then
        
    end  % channelTwo for loop  
end  % channelOne for loop

% user message
disp(['Calculated pli for ', num2str(channelNo*(channelNo-1)/2), ' channel pairings']);


return


function wpliRes = wpli(epochData)
%% Weighted Phase-Locking Index
%
% USAGE: wpliRes = wpli(epochData)
%
% Function to calculate weighted phase-locking indices (wPLIs) between a 
% set of channels / time series. IMPORTANT: Works on complex data!
%
% Input(s):
% epochData     - The input is a matrix of complexes where each row is a 
%               separate channel / time series and columns correspond to 
%               samples.  
%
% Output(s):
% pliRes        - Matrix (sized no. of channels X no. of channels) of wPLI
%               values where entry i,j is wPLI between channels i and j. As
%               the matrix is symmetric across the diagonal, there are
%               values only in the upper triangle, the rest is NaN.
%
% Relevant papers:
% Vinck et al., 2011. An improved index of phase-synchronization for 
%   electrophysiological data in the presence of volume-conduction, noise 
%   and sample-size bias. NeuroImage.
% See also:
% Stam et al., 2007. Phase lag index: Assessment of functional connectivity 
%   from multi channel EEG and MEG with diminished bias from common sources. 
%   Hum. Brain Mapp.

% IMPORTANT NOTE:
% The Vinck et al. paper defines wPLI based on cross-spectral density and 
% requires data over many trials for estimation. Here we use a version for
% a narrowband-filtered analytical signal. The idea is the same - the 
% magnitude of the "difference" (division) of complexes can be used for 
% weighting the normal PLI
%


%% Input checks

% number of input args
if nargin ~= 1
    error('Function wpli requires input arg "epochData"!');
end
% is input complex?
if isreal(epochData)
    error('Input arg "epochData" should be complex!');
end

% get number of channels and samples
channelNo = size(epochData, 1);
sampleNo = size(epochData, 2);

% user message
disp([char(10), 'Called wpli on data with ', num2str(channelNo),... 
    ' channels, each with ', num2str(sampleNo), ' samples']);


%% Loop across channel pairings

% preallocate results variable
wpliRes = nan(channelNo, channelNo);
% variable for storing channel pairings we already calculated plv for
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
            
            % we first declare the imaginary difference data we work with 
            diffData = imag(epochData(channelOne, :)./epochData(channelTwo, :));
            % then calculate wPLI
            wpliRes(channelOne, channelTwo) = abs(mean(abs(diffData).*sign(diffData)))/mean(abs(diffData));          
            
        end  % pairings if-then
        
    end  % channelTwo for loop  
end  % channelOne for loop

% user message
disp(['Calculated wpli for ', num2str(channelNo*(channelNo-1)/2), ' channel pairings']);


return


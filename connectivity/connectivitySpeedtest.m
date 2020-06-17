function [elapsedTime, connResults, pairingsNo] = connectivitySpeedtest(channelNo, epochL, measure)

%% Estimating the speed of connectivity calculations
%
% USAGE: [elapsedTime, connResults, pairingsNo] = connectivitySpeedtest(channelNo, epochL, measure)
%
% The function generates random complex data as specified in the inputs and
% then calculates PLI, PLV, iPLV or wPLI across all unique channel pairings. 
% Generated data has size of [channelNo, epochL]
%
% Input(s):
% channelNo     - number of channels to generate random phase data for
% epochL        - epoch length in samples 
% measure       - string specifying the connectivity measure:
%               'PLI', 'PLV', 'iPLV' or 'wPLI'
%
% Output(s):
% elapsedTime   - time it took for the function to finish
% connResults   - connectivity results in a matrix of size [channelNo, channelNo],
%           note that the diagonal and the lower half is filled with NaN
% pairingsNo    - number of unique channel pairings ( = no. of connectivity
%           calculations performed)
%
% Relevant papers for each connectivity measurement type:
%
% PLI:
% Stam et al., 2007. Phase lag index: Assessment of functional connectivity 
%   from multi channel EEG and MEG with diminished bias from common sources. 
%   Hum. Brain Mapp.
%
% PLV:
% Lachaux et al., 1999. Measuring phase synchrony in brain signals. 
%   Hum. Brain Mapp.
% Mormann et al., 2000. Mean phase coherence as a measure for phase 
%   synchronization and its application to the EEG of epilepsy patients. 
%   Physica D.
%
% iPLV:
% Palva, S., & Palva, J. M., 2012. Discovering oscillatory interaction 
%   networks with M/EEG: challenges and breakthroughs. Trends in cog. sci.
% Palva et al., 2018. Ghost interactions in MEG/EEG source space: A note 
%   of caution on inter-areal coupling measures. Neuroimage.
%
% wPLI:
% Vinck et al., 2011. An improved index of phase-synchronization for 
%   electrophysiological data in the presence of volume-conduction, noise 
%   and sample-size bias. NeuroImage.
%
% These papers might be of interest to you as well:
% Nolte et al., 2004. Identifying true brain interaction from EEG data 
%   using the imaginary part of coherency. Clin. Neurophys.
%
%


%% Input checks

if nargin ~= 3
    error('Need input args "channelNo", "epochL" and "measure"');
end
% sanity checking inputs
if ~ismember(channelNo, 1:1000)
    error('Are you sure about input arg "channelNo"? Prefer integer in the 1:1000 range...');
end
if ~ismember(epochL, 10:10:10^9)
    error('Are you sure about input arg "epochL"? Prefer integer in the 10:10:10^9 range...');
end
if ~ismember(measure, {'PLI', 'PLV', 'iPLV', 'wPLI'})
    error('Input arg "measure" should be one of the following: "PLI", "PLV", "iPLV" or "wPLI"');
end

% user message
disp([char(10), 'Called connectivitySpeedtest with input args:',...
    char(10), 'number of channels: ', num2str(channelNo),...
    char(10), 'epoch length in samples: ', num2str(epochL),...
    char(10), 'connectivity measure: ', measure]);


%% Basics: generate data, start timer

% generate complex data with real and imaginary parts both between -1:1 
data = (rand(epochL, channelNo)-0.5)*2*1i+(rand(epochL, channelNo)-0.5)*2;
% get phase for PLI, PLV and iPLV calculations
phaseData = angle(data);

% preallocate results matrix
connResults = nan(channelNo);

% user message
disp('Generated random data, calculating connectivity for all unique channel pairings...');

startTime = tic;


%% Loop through channel pairings, calculate PLI/PLV/iPLV/wPLI

for channelOne = 1:channelNo  
    for channelTwo = 1:channelNo
        
        % only calculate connectivity for upper triangle of connResults
        % matrix
        if channelOne < channelTwo
            
            % check connectivity measure type
            switch measure
                
                case 'PLI'
                    % PLI - relevant paper:
                    % Stam et al., 2007. Phase lag index: Assessment of functional connectivity from multi channel EEG and MEG with diminished bias from common sources. Hum. Brain Mapp.
                    connResults(channelOne, channelTwo) = abs(mean(sign(phaseData(:, channelOne)-phaseData(:, channelTwo))));
            
                case 'PLV'
                    % PLV
                    % We use the Mormann et al. version here:
                    connResults(channelOne, channelTwo) = abs(sum(exp(1i*(phaseData(:, channelOne)-phaseData(:, channelTwo))))/epochL);

                case 'iPLV'
                    % iPLV
                    % Simply the imaginary part of PLV:
                    connResults(channelOne, channelTwo) = abs(imag(sum(exp(1i*(phaseData(:, channelOne)-phaseData(:, channelTwo))))/epochL));                    
                    
                case 'wPLI'
                    % weighted PLI
                    % Note that the Vinck et al. paper defines wPLI based
                    % on cross-spectral density and requires data over many
                    % trials for estimation. Here we use a version for
                    % a narrowband-filtered analytical signal. The idea is
                    % the same - the magnitude of the "difference" 
                    % (division) of complexes can be used for weighting 
                    
                    % For readability we separately declare the imaginary
                    % difference data we work with 
                    diffData = imag(data(:, channelOne)./data(:, channelTwo));
                    
                    % wPLI
                    connResults(channelOne, channelTwo) = abs(mean(abs(diffData).*sign(diffData)))/mean(abs(diffData));      
                    
            end  % switch measure
                       
        end  % if
        
    end  % channelTwo
    
end  % channelOne


%% Finish, return 

elapsedTime = toc(startTime); 
pairingsNo = channelNo*(channelNo-1)/2;

disp(['Done!', ...
    char(10), 'Overall time elapsed: ', num2str(round(elapsedTime, 4)), ' secs',...
    char(10), 'Number of unique channel pairings: ', num2str(pairingsNo),...
    char(10), 'Elapsed time per channel pairing: ',...
    num2str(round(elapsedTime/pairingsNo, 4)), ' secs']);


return


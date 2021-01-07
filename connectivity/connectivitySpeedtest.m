function [elapsedTime, connResults, pairingsNo] = connectivitySpeedtest(metric, channelNo, sampleNo, varargin)

%% Estimating the speed of connectivity calculations
%
% USAGE: [elapsedTime, connResults, pairingsNo] = connectivitySpeedtest(metric, channelNo, sampleNo, repNo=10)
%
% The function generates random data (size as specified by inputs 
% "channelNo" and "sampleNo"), estimates connectivity ("metric") 
% across all unique channel pairings, then returns the time required for 
% connectivity estimation. 
%
% Connectivity metrics are defined outside this function. Note that all
% available metrics - as of now - are symmetric (=undirected).
%
% Mandatory inputs:
% metric       - Char array specifying the connectivity measure, one of
%               {'pli', 'wpli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}
% channelNo     - Numeric value, number of channels to generate random 
%               data for. In range 1:1000.
% sampleNo      - Numeric value, epoch length in samples. In range of
%               1:10^6.
%
% Optional input:
% repNo         - Numeric value, number of runs / repetitions for
%               estimation. In range 1:1000, defaults to 10.
%
% Outputs:
% elapsedTime   - Numeric vector, time required for connectivity
%               estimation, its values are in secs.
% connResults   - Numeric array, contains the connectivity results. Its 
%               dimensions are repNo X ROIs/channels X ROIs/channels. 
%               Note that the diagonal and the lower half is filled with 
%               NaN, only upper triangle is populated with values.
% pairingsNo    - Numeric value, number of unique channel pairings per
%               repetition ( = no. of connectivity calculations performed), 
%               equals (N-1)*N/2 where N is the number of ROIs/channels.
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
% orthAmpCorr:
% Coquelet et al., 2020. Comparing MEG and high-density EEG for intrinsic 
%   functional connectivity mapping. NeuroImage.
%
% ampCorr: 
% Standing for amplitude envelope correlation, widely used, self-evident.
%
% These papers might be of interest to you as well:
% Nolte et al., 2004. Identifying true brain interaction from EEG data 
%   using the imaginary part of coherency. Clin. Neurophys.
%


%% Input checks

if ~ismember(nargin, 3:4)
    error(['Function connectivitySpeedtest requires input args "metric", ',...
        '"channelNo" and "sampleNo" while arg "repNo" is optional!']);
end
% mandatory args
if ~ismember(metric, {'pli', 'wpli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "measure" should be one of {''pli'', ''wpli'', ''plv'', ''iplv'', ''ampCorr'', ''orthAmpCorr''}');
end
if ~ismember(channelNo, 1:1000)
    error('Are you sure about input arg "channelNo"? Prefer integer in the 1:1000 range...');
end
if ~ismember(sampleNo, 1:10^6)
    error('Inop "sampleNo"? Prefer integer in the 1:10^6 range...');
end
if ~ismember(metric, {'pli', 'wpli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "measure" should be one of {''pli'', ''wpli'', ''plv'', ''iplv'', ''ampCorr'', ''orthAmpCorr''}');
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


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
% Generated data is normally distributed pseudorandom (randn)
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
    error('Input arg "channelNo" should be a numeric value in the 1:1000 range!');
end
if ~ismember(sampleNo, 1:10^6)
    error('Input arg "sampleNo" should be a numeric value in the 1:10^6 range!');
end
% optional args
if ~isempty(varargin)
    if isnumeric(varargin{1}) && ismember(varargin{1}, 1:1000)
        repNo = varargin{1};
    end
else
    repNo = 10;
end

% user message
disp([char(10), 'Called connectivitySpeedtest with input args:',...
    char(10), 'Connectivity metric: ', metric,...
    char(10), 'Number of channels: ', num2str(channelNo),...
    char(10), 'Epoch length in samples: ', num2str(sampleNo),...
    char(10), 'Number of runs: ', num2str(repNo)]);


%% Preallocate output vars

elapsedTime = nan(repNo, 1);
connResults = nan(repNo, channelNo, channelNo);
pairingsNo = (channelNo-1)*channelNo/2;


%% Runs / repetititions loop

for repIdx = 1:repNo
    
    % Generate data
    data = randn(channelNo, sampleNo);
    
    % Get complex / phase data if the metric requires it as input
    if ismember(metric, {'pli', 'wpli', 'plv', 'iplv'})
        data = hilbert(data')';
        if ismember(metric, {'pli', 'plv', 'iplv'})
            data = angle(data);
        end
    end
    
    % timer
    startTime = tic;
    
    % connectivity estimation
    switch metric
        
        case 'pli'
            connResults(repIdx, :, :) = pli(data, 0);
        case 'wpli'
            connResults(repIdx, :, :) = wpli(data, 0);
        case 'plv'
            connResults(repIdx, :, :) = plv(data, 0);
        case 'iplv'
            connResults(repIdx, :, :) = iplv(data, 0);
        case 'ampCorr'
            connResults(repIdx, :, :) = ampCorr(data, 0);
        case 'orthAmpCorr'
            connResults(repIdx, :, :) = orthAmpCorr(data, 0);    
            
    end
    
    % get elapsed time
    elapsedTime(repIdx) = toc(startTime); 
    
end


%% User feedback, returning

disp([char(10), 'Done!', ...
    char(10), 'Overall time spent on connectivity estimation (full ', num2str(repNo), ' runs): ',... 
    num2str(round(sum(elapsedTime), 5)), ' secs',...
    char(10), 'Median time for connectivity estimation on one epoch: ',... 
    num2str(round(median(elapsedTime), 5)), ' secs',...
    char(10), 'Median time per channel pairing: ',...
    num2str(round(median(elapsedTime)/pairingsNo, 5)), ' secs']);
    
    
return


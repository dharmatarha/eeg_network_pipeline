function getMeanTimeSeries(freq, dirAngle, dirEnv, subjects)

%% Incremental averaging from polar data (angle + envelope)
%
% USAGE: getMeanTimeSeries(freq, dirAngle, dirEnv, subjects = {'s02', 's03', ...})
%
% Calculates the group average of EEG data from the polar form of 
% individual data. 
%
% Due to the amount of data to handle, the script uses incremental 
% averaging, calling only one dataset into memory and keeping a running 
% average. 
%
% Resulting group mean is saved into the current working directory and is 
% named 'FREQUENCYBAND_timeSeriesAvg.mat'.
%
% The script asssumes that the each angle and envelope .mat contains a 
% struct "EEG" with data stored in "EEG.data" (EEGlab convention, data in EEG.data). 
% Further assumption is that file naming follows the 
% 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat') convention and 
% that all relevant angle and envelope files are either in the working 
% directory or in the folders supplied in the arguments.
% 
% Mandatory inputs:
% freq       - Frequency band (string) to work with. Needs to be the same 
%        as used in the file names (e.g. 'alpha' for files like 
%        's01_alpha.mat'). One of {'alpha', 'beta', 'gamma', 'delta',
%        'theta'}
% dirAngle   - Directory path (string) pointing to the angle files. 
% dirEnv    - Directory path (string) pointing to the envelope files. 
% 
% Optional input:
% subjects   - List of subjects (cell array) whose data we average. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
%  
% NOTES: 
% (1) As with other functions in the pipeline we expect data to have 4
% dimensions: channel/ROI, sample, epoch and story/stimulus
%


%% Input checks

if ~ismember(nargin, [3, 4])
    error('Input args "freq", "dirAngle" and "dirEnv" are required!');
end
if nargin == 3
    subjects = {'s02','s03','s04','s05','s06','s07','s08','s09'...
     ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
     's21','s22','s23','s24','s25','s26','s27','s28'};
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end
if ~exist(dirAngle, 'dir')
    error('Input arg "dirAngle" is not a valid path to a folder!');
end
if ~exist(dirEnv, 'dir')
    error('Input arg "dirEnv" is not a valid path to a folder!');
end

% user message
disp([char(10), 'Starting getMeanTimeSeries function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Angle files directory: ', dirAngle,...
    char(10), 'Envelope files directory: ', dirEnv,...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% load first angle + envelope data, check dimensions
angleCheckFile = [dirAngle, '/', subjects{1}, '_', freq, '.mat'];
load(angleCheckFile, 'EEG');
angleData = EEG.data;
[roiNo, sampleNo, epochNo, stimNo] = size(angleData);
% sanity check - is it really phase data, i.e. between -pi +pi
if any(any(any(any(EEG.data>pi, 1), 2), 3), 4) || any(any(any(any(EEG.data<-pi, 1), 2), 3), 4)
    error('There is data in first angle file outside the [-pi +pi] range, is it really phase data?');
end
% check for equal sizes of angle + envelope data
envCheckFile = [dirEnv, '/', subjects{1}, '_', freq, '.mat'];
load(envCheckFile, 'EEG');
envData = EEG.data;
if ~isequal(size(envData), [roiNo, sampleNo, epochNo, stimNo])
    error('Envelope data has different size than corresponding angle data for the first files!');
end

% user message
disp([char(10), 'First data files had ', num2str(roiNo), ' channels/ROIs, ',... 
    num2str(stimNo), ' different stimuli (stories), ',...
    num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
    ' samples for each epoch. Assuming same dimensions for each subsequent data file.']);


%% Loop to calculate means

disp([char(10), 'Starting calculation of group mean...', char(10)])

% for creating the running average, we simply start with the first
% subject's data loaded in the previous block
[realData, ~] = pol2cart(angleData, envData);  % transform polar to cartesian and get rid of the hilbert part
runningAvg = realData;

% user message
disp(['Finished with subject ', subjects{1}]);

% we start the subject loop from the second one
for subIndex = 2:subNo
    
    % load next angle + envelope data
    nextAngleFile = [dirName, '/', subjects{subIndex}, '_', freq,'.mat'];
    load(nextAngleFile);
    angleData = EEG.data;
    nextEnvFile = [dirName, '/', subjects{subIndex}, '_', freq,'.mat'];
    load(nextEnvFile);
    envData = EEG.data;

    % sanity check - data sizes
    if ~isequal(size(angleData), [roiNo, sampleNo, epochNo, stimNo]) || ~isequal(size(envData), [roiNo, sampleNo, epochNo, stimNo])
        error(['Angle or envelope data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
    end       
    % sanity check - is it really phase data, i.e. between -pi +pi
    if any(any(any(any(angleData>pi, 1), 2), 3), 4) || any(any(any(any(angleData<-pi, 1), 2), 3), 4)
        error(['Angle data for subject ', subjects{subIndex}, ' has values outside the [-pi +pi] range, is it really phase data?']);
    end
    
    % transform to real time series
    [realData, ~] = pol2cart(angleData, envData);
      
    % incremental averaging
    runningAvg = runningAvg + (realData - runningAvg)/subIndex;

    % user message
    disp(['Finished with subject ', subjects{subIndex}]);    
    
end


%% save out the results

saveFile = [pwd, '/', freq, '_timeSeriesAvg.mat'];
save(saveFile, 'runningAvg', 'subjects');

disp([char(10), 'Finished with real time series group mean', ... 
    ', saved out results to ', saveFile, char(10)]);


return












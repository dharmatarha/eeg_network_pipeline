function selectConnEpochs(dirName, method, freq)
%% Function for randomly selecting epochs in the EEG resting state dataset
%
% USAGE: selectConnEpochs(dirName, method, freq)
% 
% The EEG resting state dataset - after preprocessing - consists of 4 sec
% epochs with 2 secs overlap, with a variable number of epochs per
% participant. Connectivity is calculated for these overlapping epochs,
% resulting in one channel/ROI no. X channel/ROI no. matrix per epoch for each
% subject.
%
% This function works on the connectivity data. It first selects 
% non-overlapping epochs (= every second) and then randomly selects 
% "N" epochs from each subject where "N" is smallest number of epochs 
% across all participants.
%
% Selected epochs for each subject are arranged into a 4D array with
% dimensions subjects X epochs X channels/ROIs X channels/ROIs. The 4D
% array is saved out into a .mat file in "dirName", called
% "group_FREQUENCYBAND_METHOD.mat".
%
% All input files are assumed to conform to the naming convention
% "SUBJECTID_FREQUENCYBAND_METHOD.mat" where method refers to the
% connectivity measure used. Files should contain a var "connRes" with the
% connectivity data sized epoch no. X channels/ROIs X channels/ROIs.
%
% Inputs:
% dirName   - Char array, path to folder containing files 
% method    - Char array, one of {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}.
% freq      - Char array, one of 
%           {'delta', 'theta', 'alpha', 'beta', 'gamma'}
%
% Output:
% group_FREQUENCYBAND_METHOD.mat saved into "dirName", containing: 
% connData  - 4D numeric array with dimensions subjects X epochs X 
%           channels/ROIs X channels/ROIs
% subjects  - Cell array, contains a subject IDs corresponding to the first
%           dimension of "connData".
% epochIndices - Cell array, contains the indices of the epochs selected
%           from the connectivity data. 
%           NOTE: Indices from the list of every second epochs!
% 


%% Input checks

if nargin ~= 3 
    error('Function selectConnEpochs requires input args "dirName", "method" and "freq"!');
end
if ~exist(dirName, 'dir')
    error('Input arg "dirName" is not a valid folder path!');
end
if ~ismember(method, {'iplv', 'plv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "method" has an unexpected value!');
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end
% standardize dir paths regarding trailing '/'
if dirName(end) == '/'
    dirName = dirName(1:end-1);
end


%% List files

% list everything matching the supplied frequency band & connectivity 
% method in the target dir
subFiles = dir([dirName, '/*_', freq, '_', method, '.mat']); 

% number of subjects
subNo = length(subFiles);

% var storing subject ids (file name part before "_freq")
subIDs = cell(subNo, 1);

% create file paths from structs
subFilePaths = cell(subNo, 1);
for i = 1:subNo
    subFilePaths{i} = [subFiles(i).folder, '/', subFiles(i).name];
    % subject ID is the part before _FREQUENCYBAND
    subIDs{i} = subFiles{i}.name(1:strfind(subFiles{i}.name, freq)-2);
end

% user message
disp([char(10), 'Found ', num2str(subNo),... 
    ' .mat files matching the supplied frequency band and ',...
    'conenctivity method in ', dirName]);
disp('Subject IDs for the files:');
disp(subIDs);


%% Load files, find the smallest epoch number

% var to store epoch numbers for each subject
epochNumbers = nan(subNo, 1);
% struct to store all loaded data - we do not check for excessive memory
% usage, that is the user's responsibility!
subData = struct('connRes', nan(100, 100, 100));
subData(subNo).connRes = nan(100, 100, 100);  % preallocate with last element of struct

% determine channel/ROI numbers in first file, use that for sanity checks
% in case of all other files
tmp = load(subFilePaths{1});
[~, refRoiNo1, refRoiNo2] = size(tmp.connRes);
% sanity check
if refRoiNo1 ~= refRoiNo2
    error(['2nd and 3rd dims of first data file at ',... 
        subFilePaths{1}, ' are not equal!']);
end

% loop through files, load them
for i = 1:subNo
    % load file content
    tmp = load(subFilePaths{i});
    % only store every second epoch!
    subData(i).connRes = tmp.connRes(1:2:end, :, :);
    % get size
    [subEpochNo, subRoiNo1, subRoiNo2] = size(subData(i).connRes);
    % sanity check
    if subRoiNo1 ~= subRoiNo2 || subRoiNo1 ~= refRoiNo1
        error(['2nd and 3rd dims of data file at ',... 
        subFilePaths{i}, ' have unexpected size, investigate!']);
    end
    epochNumbers(i) = subEpochNo;
end
    

%% Determine cutoff 

% report epoch numbers
disp([char(10), 'Epoch numbers per subjects in ascending order: ']);
disp(sort(epochNumbers));

% ask for input for cutoff point
question = [char(10), 'What should be the number of epochs kept (cutoff) ',...
    'for subjects? \n If your choice is larger than the minimal epoch number, ',...
    'subjects with less epochs than the cutoff will be left out of the final results array!', char(10)];
cutoffStr = input(question, 's');

% get numeric, check for number of remaining subjects
cutoff = str2double(cutoffStr);
if cutoff > min(epochNumbers)
    subLeft = sum(epochNumbers >= cutoff);
    disp([char(10), 'Cutoff is ', num2str(cutoff),... 
        '. This is larger than the minimal epoch number, ',... 
        num2str(subLeft), ' out of ', num2str(subNo),... 
        ' subjects will remain in final data array.']);
else
    subLeft = subNo;
    disp([char(10), 'Cutoff is equal / smaller than the minimal ',...
        'epoch number, all subjects will remain in final data array.']);
end


%% Loop through data, select epochs

% preallocate
connData = zeros(subLeft, cutoff, refRoiNo1, refRoiNo1);
subjects = cell(subLeft, 1);
epochIndices = cell(subLeft, 1);

% loop through data sets
subCounter = 0;
for i = 1:subNo

    % survives cutoff?
    if epochNumbers(i) >= cutoff
        subCounter = subCounter + 1;
        
        % collect subject ID
        subjects{subCounter} = subIDs{i};
        
        % need to sample from epochs?
        if epochNumbers(i) > cutoff
            
            % random epoch indices we use for epoch selection
            epochIndices{subCounter} = randperm(epochNumbers(i), cutoff);
            % selected epochs into output array
            connData(subCounter, :, :, :) = subData(i).connRes(epochIndices{subCounter}, :, :);
        
        % if cutoff equals the number of epochs, just copy connectivity
        % data
        elseif epochNumbers(i) == cutoff
            
            connData(subCounter, :, :, :) = subData(i).connRes;
            
        end  % if epochNumbers > cutoff
        
    end  % if epochNumbers >= cutoff
    
end  % for subNo
    
% user message    
disp([char(10), 'Selected epochs from all eligible subjects (',...
    num2str(subCounter), ' participants).']);


%% Save output vars

saveF = [dirName, '/group_', freq, '_', method, '.mat'];
save(saveF, 'connData', 'subjects', 'epochIndices');

% user message
disp([char(10), 'Saved out results, returning...']);


return


















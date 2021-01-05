function selectConnEpochs(dirName, method, freq, varargin)
%% Function for randomly selecting epochs in the EEG resting state dataset
%
% USAGE: selectConnEpochs(dirName, method, freq, subjects=[], epochIndices=[])
% 
% Function for selecting an equal number of epochs from a set of subjects
% with variable epoch numbers.
%
% The EEG resting state dataset - after preprocessing - consists of 4 sec
% epochs with 2 secs overlap, with a variable number of epochs per
% participant. Connectivity is calculated for these overlapping epochs,
% resulting in a matrix sized [ROI no. X ROI no.] per epoch for each
% subject.
%
% This function works on the connectivity matrices. It first selects 
% non-overlapping epochs (= every second) and then randomly selects 
% "N" epochs from each subject where "N" is the smallest number of epochs 
% across all participants. "N" is determined via user input after
% displaying the epoch numbers detected.
%
% NOTE THAT EPOCH INDICES REFER TO POSITION AFTER SELECTING EVERY SECOND
% EPOCH!
%
% Selected epochs for each subject are arranged into a 4D array with
% dimensions subjects X epochs X ROIs X ROIs. The 4D
% array is saved out into a .mat file in "dirName", called
% "group_FREQUENCYBAND_METHOD.mat".
%
% All input files are assumed to conform to the naming convention
% "SUBJECTID_FREQUENCYBAND_METHOD.mat" where method refers to the
% connectivity measure used. Files should contain a var "connRes" with the
% connectivity data sized [epoch no. X ROI no. X ROI no.].
%
% Mandatory inputs:
% dirName   - Char array, path to folder containing files 
% method    - Char array, one of {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}.
% freq      - Char array, one of 
%           {'delta', 'theta', 'alpha', 'beta', 'gamma'}
%
% Optiional inputs:
% subjects  - Cell array, includes the subject IDs for the set of subjects
%           to include (char array in each cell). Determines the sequence 
%           of subjects in output vars "connData" and "epochIndices" 
%           (unless "epochIndices" is also supplied as input arg). 
%           Defaults to [], meaning that all subjects with connectivity 
%           data files in "dirName" matching "method" and "freq" are taken 
%           into account.
% epochIndices - Cell array sized [length(subjects, 1]. Each cell contains
%           a numeric vector corresponding to the epoch indices to select 
%           from the corresponding subject's data. Defaults to [], meaning
%           that epochs are randomly selected for each subject. If
%           supplied, its length must match that of "subjects".
%
% Output:
% group_FREQUENCYBAND_METHOD.mat saved into "dirName", containing: 
% connData  - 4D numeric array with dimensions subjects X epochs X 
%           ROIs X ROIs. Contains connectivity data from selected epochs.
% subjects  - Cell array, contains subject IDs corresponding to the first
%           dimension of "connData".
% epochIndices - Cell array, contains the indices of the epochs selected
%           from the connectivity data. Empty if all epochs were included.
%           NOTE: Indices from the list of every second epochs!
% 


%% Input checks

if ~ismember(nargin, 3:5) 
    error(['Function selectConnEpochs requires input args "dirName", ',...
        '"method" and "freq" while args "subjects" and "epochIndices" are optional!']);
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
% check optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if iscell(varargin{v}) && all(cellfun(@ischar, varargin{v})) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif iscell(varargin{v}) && all(cellfun(@isnumeric, varargin{v})) && all(cellfun(@isvector, varargin{v})) && ~exist('epochIndices', 'var')
            epochIndices = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to optional args "subjects" or "epochIndices"!');
        end
    end
end
% defaults
if ~exist('subjects', 'var')
    subjects = [];
end
if ~exist('epochIndices', 'var')
    epochIndices = [];
end
% extra checks
if ~isequal(length(subjects), length(epochIndices))
    error('Input args "subjects" and "epochIndices" must have equal length!');
end
if ~isempty(epochIndices)
    epochNo1 = length(epochIndices{1});
    epochIdxLengths = cellfun('length', epochIndices);
    if ~all(ismember(epochIdxLengths, epochNo1))
        error('Not all cells in "epochIndices" contain the same number of indices (=differently sized numeric arrays in the cells)!');
    end
end
% standardize dir paths regarding trailing '/'
if dirName(end) == '/'
    dirName = dirName(1:end-1);
end

% user message
disp([char(10), 'Called function selectConnEpochs with input args: ',...
    char(10), 'Directory for connectivity data files: ', dirName, ...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Supplied subject list: ', num2str(~isempty(subjects)),...
    char(10), 'Supplied epoch indices: ', num2str(~isempty(epochIndices))]);


%% List files

% if there was no subject list supplied as input arg, check "dirName" for
% files matching "freq" and "method"
if isempty(subjects)
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
        subIDs{i} = subFiles(i).name(1:strfind(subFiles(i).name, freq)-2);
    end

    % user message
    disp([char(10), 'Found ', num2str(subNo),... 
        ' .mat files matching the supplied frequency band, metric ',... 
        char(10), 'and connectivity metric in ', dirName]);
    disp('Subject IDs for the files:');
    disp(subIDs);

% if there was a subject list supplied, check if the corresponding files exist    
else
    subNo = length(subjects);
    subFilePaths = cell(subNo, 1);  % preallocate
    % create file paths corresponding to supplied subject IDs, frequency
    % band and connectivity metric
    for i = 1:subNo
        subFilePaths{i} = [dirName, '/', subjects{i}, '_', freq, '_', method, '.mat'];
    end
    % check if these files exist
    pathValidity = cellfun(@exist, subFilePaths);
    % if there is any file missing, print the list of missing files and
    % error out
    if ~any(pathValidity)
        warning('Could not find connectivity data for all subjects, these files are missing: ');
        disp(subFilePaths(~pathValidity));
        error('Cannot continue without a valid file list for connectivity data...');
    end
    
    % user message
    disp([char(10), 'Found ', num2str(subNo),... 
        ' .mat files matching the supplied subject IDs, ',...
        char(10), 'frequency band and connectivity metric in ', dirName]);
end


%% Load files, get epoch numbers

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % !!! IMPORTANT HARDCODED METHOD !!!
    % only store every second epoch!
    subData(i).connRes = tmp.connRes(1:2:end, :, :);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get size
    [subEpochNo, subRoiNo1, subRoiNo2] = size(subData(i).connRes);
    % sanity check
    if subRoiNo1 ~= subRoiNo2 || subRoiNo1 ~= refRoiNo1
        error(['2nd and 3rd dims of data file at ',... 
        subFilePaths{i}, ' have unexpected size, investigate!']);
    end
    epochNumbers(i) = subEpochNo;
end
    
    
%% Determine cutoff for epoch numbers to include for each subject - only if epoch indices were not supplied

if isempty(epochIndices)
    
    %% Determine cutoff 

    % report epoch numbers
    disp([char(10), 'Epoch numbers per subjects in descending order: ']);
    disp(sort(epochNumbers, 'descend'));

    % ask for input for cutoff point
    question = [char(10), 'What should be the number of epochs kept (cutoff) for subjects? \n',... 
        'If your choice is larger than the minimal epoch number, \n',...
        'subjects with less epochs than the cutoff \n',... 
        'will be left out of the final results array!', char(10)];
    cutoffStr = input(question, 's');

    % get numeric, check for number of remaining subjects
    cutoff = str2double(cutoffStr);
    if cutoff > min(epochNumbers)
        subLeftNo = sum(epochNumbers >= cutoff);
        disp([char(10), 'Cutoff is ', num2str(cutoff), '. '... 
            char(10), 'This is larger than the minimal epoch number, ',... 
            char(10), num2str(subLeftNo), ' out of ', num2str(subNo),... 
            ' subjects will remain in final data array.']);
    else
        subLeftNo = subNo;
        disp([char(10), 'Cutoff is equal / smaller than the minimal epoch number, ',... 
            char(10), 'all subjects will remain in final data array.']);
    end
   
% if there were epoch indices supplied, check if they are compatible 
% with epoch numbers detected, error out if not
else
    indexMatch = epochIdxLengths <= epochNumbers;
    if ~all(indexMatch)
        warning('More indices in "epochIndices" than available epochs for entries: ');
        disp(find(~indexMatch));
        error('The number of indices in "epochIndices" is larger than the number of available epochs for at least one subject!')
    end
end


%% Loop through data, select epochs

% preallocate result vars
if isempty(subjects)
    subLeftID = cell(subLeftNo, 1);
end
if isempty(epochIndices)
    selectedEpochs = cell(subLeftNo, 1);
    connData = zeros(subLeftNo, cutoff, refRoiNo1, refRoiNo1);
else
    connData = zeros(subNo, epochNo1, refRoiNo1, refRoiNo1);
end

% loop through data sets
subCounter = 0;
for i = 1:subNo
    
    % if there were no epoch indices supplied, we randomly select the
    % "cutoff" number of epochs for each subject with enough epochs
    if isempty(epochIndices)
        
        % subject survives cutoff?
        if epochNumbers(i) >= cutoff
            subCounter = subCounter + 1;

            % collect subject ID
            subLeftID{subCounter} = subIDs{i};

            % need to sample from epochs?
            if epochNumbers(i) >= cutoff

                % random epoch indices we use for epoch selection
                selectedEpochs{subCounter} = randperm(epochNumbers(i), cutoff);
                % selected epochs into output array
                connData(subCounter, :, :, :) = subData(i).connRes(selectedEpochs{subCounter}, :, :);

            end  % if epochNumbers > cutoff

        end  % if epochNumbers >= cutoff
    
    % if there were epoch indices supplied, just select those epochs 
    % and copy them into output array    
    else
        connData(i, :, :, :) = subData(i).connRes(epochIndices{i}, :, :);
        
    end  % if isempty(epochIndices)
           
end  % for subNo

% epochIndices and subjects are alos output vars
if isempty(subjects)
    subjects = subLeftID;
end
if isempty(epochIndices)
    epochIndices = selectedEpochs;
end

% user message    
disp([char(10), 'Selected epochs from all eligible subjects (',...
    num2str(length(subjects)), ' participants).']);


%% Save output vars

saveF = [dirName, '/group_', freq, '_', method, '.mat'];
save(saveF, 'connData', 'subjects', 'epochIndices');

% user message
disp([char(10), 'Saved out results, returning...']);


return













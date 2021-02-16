function [subjectFiles, selectedEpochs, epochData] = selectEpochs(dirName, filePattern, varName, varargin)
%% Function for randomly selecting epochs in multi-subject EEG dataset
%
% USAGE: [subjectFiles, selectedEpochs, epochData] = selectEpochs(dirName, 
%                                                               filePattern, 
%                                                               varName, 
%                                                               epochNo=[], 
%                                                               subjects=[], 
%                                                               epochMask=[], 
%                                                               epochIndices=[])
% 
% Function for selecting (randomly) an equal number of epochs from a set 
% of subjects with variable epoch numbers.
%
% Basic usage is pointing the function to a set of files using input
% args "dirName", "filePattern" and "varName", and then using the 
% command line dialogue for selecting the number of epochs. 
%
% Each file matching "filePattern" will be treated as one subject's data.
% Data is expected to be in a 3D numeric array named "varName", with the 
% first dimension marking epochs (e.g. sized [epochs X ROIs X ROIs] for 
% connectivity data). 
%
% If "epochNo" is specified, "epochNo" number of epochs are selected from
% each subject's data. If a subject has less than "epochNo" epochs, her
% data is left out from the final data set.
%
% If "subjects" is specified, we try to match the char arrays in "subjects"
% to the beginning of the file names and only include the files we have 
% matching "subject" cells for.
%
% If "epochMask" is specified, we only take into account the epochs
% in the mask (values "true" or "1" in the mask mark included epochs).
% Arg "epochMask" might be longer or shorter than the number of epochs for 
% one or more subjects, we still apply it. Epochs outside the length of the
% mask are not considered.
%
% If "epochIndices" are specified, its cells contain the epochs to be
% selected for each subject. Useful for consistent selection of epochs
% across different aspects of the data (e.g. different frequency bands). 
%
% Selected epochs for each subject are arranged into a 4D array with
% dimensions subjects X epochs X ROIs X ROIs (output "epochData").
% 
% IMPORTANT! The script first loads the files one-by-one and collects the
% number of available epochs for each subject. This part might take a few
% minutes for larger data sets.
%
% Mandatory inputs:
% dirName           - Char array, path to folder containing files 
% filePattern       - Char array, file name part for specifying the files 
%                   to include in the selection process. The function uses
%                   "dir" for finding files, so asterisk wildcards are
%                   allowed. An extra ".mat" ending is always added when
%                   listing the files, do not include that in "filePattern"
% varName           - Char array, name of the variable containing the data
%                   of interest
%
% Optional inputs:
% epochNo           - Numeric value, specifies the number of epochs we
%                   select randomly for each subject. If there are less
%                   than "epochNo" number of epochs for a subject, that
%                   subject's data is left out (no epochs are selected). If
%                   not specified, the function informs the user about the
%                   the epoch numbers found in the files and then prompts 
%                   the user for an "epochNo" value (a cutoff) to use.
% subjects          - Cell array of char arrays, includes the subject IDs 
%                   for the set of subjects to include. File names are 
%                   assumed to start with the subject ID (e.g. 
%                   'l3_s10_alpha_plv.mat' for subject ID 'l3_s10'). 
%                   Determines the sequence of subjects in the output array 
%                   "eopchData.
%                   If not specified, the function takes all files matching
%                   the mandatory input args into account. 
% epochMask         - Numeric (binary) or logical vector. Specifies which
%                   epochs to take into account for random selection. E.g.
%                   if epochMask = [0 1 1 0 1], only the 2nd, 3rd and 5th
%                   epochs are used for selection, irrespective of how many
%                   epochs are actually available for any subject. Useful
%                   e.g. for selecting non-overlappping epochs if epochs
%                   are otherwise overlapping. If not specified, all epochs
%                   are taken into account.
% epochIndices      - Cell array of numeric vectors. Its length must match
%                   the number of subjects / data files taken into account.
%                   Each cell contains a numeric vector corresponding to 
%                   the epoch indices to select from the corresponding 
%                   subject's data. If not specified, epochs are randomly 
%                   selected for each subject.
%
% Outputs:
% subjectFiles      - Cell array of char arrays. Contains the names of all
%                   files taken into account for epoch selection. Its
%                   length equals the size of the first dimension of
%                   "epochData" and the length of "selectedEpochs".
% selectedEpochs    - Cell array of numeric vectors. Each cell contains a
%                   numeric vector corresponding to the epoch indices of
%                   selected epochs. E.g. if selectedEpochs{2} is [3 10
%                   20], it means that for the second subject
%                   (subjectFiles{2}) "epochData" contains the 3rd, 10th
%                   and 20th epochs. Its length equals the length of
%                   "subjectFiles" and the first dimension of "epochData".
%                   The length of each numeric vector is the number of
%                   selected epochs.
% epochData         - Numeric array, 4D, firs dimension is subjects, second 
%                   is epochs, the rest is the EEG data of the selected 
%                   epochs sized. E.g. for connectivity data it is sized 
%                   [subjects X epochs X ROIs X ROIs].
%
%
%
% Originally intended use case:
% The EEG resting state dataset - after preprocessing - consists of 4 sec
% epochs with 2 secs overlap, with a variable number of epochs per
% participant. Connectivity is calculated for these overlapping epochs,
% resulting in a matrix sized [ROI no. X ROI no.] per epoch for each
% subject. Connectivity data is stored in files with a naming convention
% SUBJECTID_FREQUENCYBAND_METHOD.mat.
% For selecting an equal number of non-overlappipng epochs for all 
% subjects, using the alpha frequency band and "plv" as the connectivity 
% estimate, we do:
% 
% % get subject ids for the resting state dataset from /utils
% load('restingStateSubjects.mat', 'subjectsRS');
% % define the mask for non-overlapping epochs (every second in our case)
% epochMask = zeros(1,300); epochMask(1:2:end) = 1;  % 300 is definitely more than the max number of epochs we have
% % call selecteEpochs
% [subjectFiles, selectedEpochs, epochData] = selectEpochs('alpha/', 'alpha_plv', subjectsRS, epochMask);
% 
%


%% Input checks

% no. if inputs
if ~ismember(nargin, 3:7) 
    error(['Function selectConnEpochs requires input args "dirName", ',...
        '"filePattern" and "varName" while args "epochNo", "subjects", ',...
        '"epochMask" and "epochIndices" are optional!']);
end
% mandatory args
if ~exist(dirName, 'dir')
    error('Input arg "dirName" is not a valid folder path!');
end
if ~ischar(filePattern)
    error('Input arg "filePattern" should be char array!');
end
if ~ischar(varName)
    error('Input arg "varName" should be char array!');
end
% optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('epochNo', 'var')
            epochNo = varargin{v}; 
        elseif iscell(varargin{v}) && all(cellfun(@ischar, varargin{v})) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif (isnumeric(varargin{v}) || islogical(varargin{v})) && isvector(varargin{v}) && numel(varargin{v})>1 && ~exist('epochMask', 'var')
            epochMask = varargin{v};
        elseif iscell(varargin{v}) && all(cellfun(@isnumeric, varargin{v})) && all(cellfun(@isvector, varargin{v})) && ~exist('epochIndices', 'var')
            epochIndices = varargin{v};
        else
            error(['At least one input arg could not be mapped nicely to ',...
                'optional args "epochNo", "subjects", "epochMask" or "epochIndices"!']);
        end
    end
end
% epochMask: if numeric, check if if also binary
if exist('epochMask', 'var') && isnumeric(epochMask) && ~all(ismember(epochMask, [0 1]))
    error('Input arg "epochMask" is numeric and not binary! Make sure it only contains zeros and ones!');
end
% defaults
if ~exist('epochNo', 'var')
    epochNo = [];
end
if ~exist('subjects', 'var')
    subjects = [];
end
if ~exist('epochMask', 'var')
    epochMask = [];
end
if ~exist('epochIndices', 'var')
    epochIndices = [];
end
% extra checks
% if "epochNo" and "epochMask" are specified, check no. of epochs in mask
if ~isempty(epochNo) && ~isempty(epochMask)  && sum(epochMask) < epochNo
    error(['Incompatible args "epochNo" and "epochMask" - number of ',...
        'epochs for selection ("epochNo") is larger than the number ',...
        'of epochs in the mask ("epochMask")!']);
end
% if "epochIndices" is specified, error out if "epochMask" or
% "epochNo" are also specified
if ~isempty(epochIndices) && (~isempty(epochNo) || ~isempty(epochMask))
    error(['Input arg "epochIndices" is specified, but "epochNo" or "epochMask" ',...
        'is also supplied. This makes no sense as "epochIndices" determines ',...
        'the epochs for selection completely.']);
end
% if "subjects" and "epochIndices" are specified, check length of cell
% arrays
if ~isempty(subjects) && ~isempty(epochIndices) && ~isequal(length(subjects), length(epochIndices))
    error('Input args "subjects" and "epochIndices" must have equal length!');
end
% if "epochIndices" contains cells with different numbers of epoch indices
if ~isempty(epochIndices)
    tmp = length(epochIndices{1});
    epochIdxLengths = cellfun('length', epochIndices);
    if ~all(ismember(epochIdxLengths, tmp))
        error('Not all cells in "epochIndices" contain the same number of indices (=differently sized numeric arrays in the cells)!');
    end
end

% standardize dir paths regarding trailing '/'
if dirName(end) == '/'
    dirName = dirName(1:end-1);
end

% user message
disp([char(10), 'Called function selectConnEpochs with input args: ',...
    char(10), 'Folder of data files: ', dirName, ...
    char(10), 'File name part for file selection: ', filePattern,...
    char(10), 'Variable name containing data: ', varName,...
    char(10), 'Supplied number of epochs to select: ', num2str(~isempty(epochNo)),...
    char(10), 'Supplied subject list: ', num2str(~isempty(subjects)),...
    char(10), 'Supplied epoch mask: ', num2str(~isempty(epochMask)),...    
    char(10), 'Supplied epoch indices: ', num2str(~isempty(epochIndices))]);


%% List files

% Check "dirName" and "filePattern" for data files
subFiles = dir([dirName, '/*', filePattern, '*.mat']); 
% number of subjects
subNo = length(subFiles);
% var storing subject file names
subjectFiles = extractfield(subFiles, 'name');

% user message
disp([char(10), 'Found ', num2str(subNo),... 
    ' .mat files matching "dirName" and "filePattern"']);

% if there was a subject list supplied, check if the corresponding files
% exist, and there is only one file per cell in "subjects"
if ~isempty(subjects)
    tmpSum = zeros(size(subjectFiles));  % preallocate
    % go through subject ids, check if there is one file for each subject
    % id
    for s = 1:subNo
        tmp = startsWith(subjectFiles, subjects{s});
        if sum(tmp)~=1
            error(['Error at matching data files to supplied subject ids. ',...
                'There was either multiple files or no file starting with ', subjects{s}]);
        end
        tmpSum = tmpSum+tmp;  % accumulate matches
    end
    % check if there was any file with multiple matches
    if any(tmpSum>1)
        disp(subjectFiles(tmpSum>1)');
        error(['Error at matching data files to supplied subject ids. ',...
            'At least one data file had more than one match to subject ids. ',...
            'See list of affected files above.']);
    end
    % if everything was fine, delete data files without a subject id match
    subjectFiles(tmpSum==0) = [];
    
    % user message
    disp([char(10), 'After matching potential data files to subject ids, ',...
        'there are ', num2str(length(subjectFiles)), ' files remaining.']);
    
end

% create file paths from structs
subFilePaths = cell(subNo, 1);
for i = 1:subNo
    subFilePaths{i} = [dirName, '/', subjectFiles{i}];
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
    connData = zeros(subNo, tmp, refRoiNo1, refRoiNo1);
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













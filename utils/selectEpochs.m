function [selectedFiles, selectedEpochs, epochData] = selectEpochs(dirName, filePattern, varName, varargin)
%% Function for randomly selecting epochs in multi-subject EEG dataset
%
% USAGE: [selectedFiles, selectedEpochs, epochData] = selectEpochs(dirName, 
%                                                               filePattern, 
%                                                               varName, 
%                                                               epochDim=1,
%                                                               epochNo=[], 
%                                                               subjects=[], 
%                                                               epochMask=[], 
%                                                               epochIndices=[],
%                                                               loadBehav='incremental')
% 
% Function for selecting (randomly) an equal number of epochs from a set 
% of subjects with variable epoch numbers.
%
% Basic usage is pointing the function to a set of files using input
% args "dirName", "filePattern" and "varName", and then using the 
% command line dialogue for selecting the number of epochs. 
%
% Each file matching "filePattern" will be treated as one subject's data.
% By default, data is expected to be in a 3D numeric array named "varName", 
% with the first dimension marking epochs (e.g. sized 
% [epochs X ROIs X ROIs] for connectivity data). If epochs are marked by
% another dimension (not first), you can specify that with "epochDim".
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
% Arg "loadBehav" specifies if all data should be loaded in one sweep
% ('onesweep') and kept in memory for epoch selection, or if loading and
% epoch selection should be done incrementally ('incremental'). The later
% is useful if a relatively small subset of epochs is selected from a large
% dataset. 
% IMPORTANT! If "loadBehav" is set to 'incremental' and neither "epochNo" 
% or "epochIndices" are specified, the files are loaded twice, one-by-one 
% in both cases - first to determine epoch numbers and then to collect 
% selected epochs. This can take a while for a large dataset.
%
% Selected epochs for each subject are arranged into a 4D array with
% dimensions subjects X epochs X ROIs X ROIs (output "epochData").
% 
% Optional args are inferred from types and values. If that is ambigious,
% position is taken into account (only in the case of differentiating
% "epochDim" and "epochNo").
%
%
% Mandatory inputs:
% dirName           - Char array, path to folder containing files 
% filePattern       - Char array, file name part for specifying the files 
%                   to include in the selection process. The function uses
%                   "dir" for finding files, so asterisk wildcards are
%                   allowed.
% varName           - Char array, name of the variable containing the data
%                   of interest
%
% Optional inputs:
% epochDim          - Numeric value, specifies which dimension of the data 
%                   array in "varName" marks epochs. Defaults to 1. 
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
% loadBehav         - Char array, one of {'onesweep', 'incremental'}.
%                   Determines if data files are loaded all at once and 
%                   kept in memory for epoch selection ('onesweep') or if
%                   only one subject's data is kept in memory at any time
%                   for epoch selection ('incremental'). The former is 
%                   memory-hungry while the latter can take longer.
%                   Defaults to 'incremental' to avoid out-of-memmory
%                   problems.
%
% Outputs:
% selectedFiles     - Cell array of char arrays. Contains the names of all
%                   files taken into account for epoch selection. Its
%                   length equals the size of the first dimension of
%                   "epochData" and the length of "selectedEpochs".
% selectedEpochs    - Cell array of numeric vectors. Each cell contains a
%                   numeric vector corresponding to the epoch indices of
%                   selected epochs. E.g. if selectedEpochs{2} is [3 10
%                   20], it means that for the second subject
%                   (selectedFiles{2}) "epochData" contains the 3rd, 10th
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
% % call selectEpochs
% [subjectFiles, selectedEpochs, epochData] = selectEpochs('alpha/', 'alpha_plv.mat', subjectsRS, epochMask);
% 
%


%% Input checks

% no. if inputs
if ~ismember(nargin, 3:9) 
    error(['Function selectConnEpochs requires input args "dirName", ',...
        '"filePattern" and "varName" while args "epochDim", "epochNo", "subjects", ',...
        '"epochMask", "epochIndices" and "loadBehav" are optional!']);
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
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('epochDim', 'var')
            epochDim = varargin{v}; 
        elseif isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('epochNo', 'var')
            epochNo = varargin{v};            
        elseif iscell(varargin{v}) && all(cellfun(@ischar, varargin{v})) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif (isnumeric(varargin{v}) || islogical(varargin{v})) && isvector(varargin{v}) && numel(varargin{v})>1 && ~exist('epochMask', 'var')
            epochMask = varargin{v};
        elseif iscell(varargin{v}) && all(cellfun(@isnumeric, varargin{v})) && all(cellfun(@isvector, varargin{v})) && ~exist('epochIndices', 'var')
            epochIndices = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'onesweep', 'incremental'}) && ~exist('loadBehav', 'var')
            loadBehav = varargin{v};             
        else
            error(['At least one input arg could not be mapped nicely to ',...
                'optional args "epochDim", "epochNo", "subjects", "epochMask", "epochIndices" or "loadBehav"!']);
        end
    end
end
% epochMask: if numeric, check if also binary
if exist('epochMask', 'var') && isnumeric(epochMask) && ~all(ismember(epochMask, [0 1]))
    error('Input arg "epochMask" is numeric and not binary! Make sure it only contains zeros and ones!');
end
% defaults
if ~exist('epochDim', 'var')
    epochDim = 1;
end
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
if ~exist('loadBehav', 'var')
    loadBehav = 'incremental';
end
% extra checks
% if "epochNo" and "epochMask" are specified, check no. of epochs in mask
if ~isempty(epochNo) && ~isempty(epochMask) 
    if sum(epochMask) < epochNo
        error(['Incompatible args "epochNo" and "epochMask" - number of ',...
            'epochs for selection ("epochNo") is larger than the number ',...
            'of epochs in the mask ("epochMask")!']);
    end
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

% make sure epochMask is a column vector
if ~isempty(epochMask) && isrow(epochMask)
    epochMask = epochMask';
end

% user message
disp([char(10), 'Called function selectConnEpochs with input args: ',...
    char(10), 'Folder of data files: ', dirName, ...
    char(10), 'File name part for file selection: ', filePattern,...
    char(10), 'Name of variable containing data of interest: ', varName,...
    char(10), 'Supplied number of epochs to select: ', num2str(~isempty(epochNo)),...
    char(10), 'Supplied subject list: ', num2str(~isempty(subjects)),...
    char(10), 'Supplied epoch mask: ', num2str(~isempty(epochMask)),...    
    char(10), 'Supplied epoch indices: ', num2str(~isempty(epochIndices)),...
    char(10), 'Data loading behavior: ', loadBehav]);


%% List files (with their paths) matching the supplied input args

[filePaths, ~] = listMatchingFiles(dirName, filePattern, subjects);


%% Load files, get epoch numbers

% user message
disp([char(10), 'Loading all files, collecting number of epochs from each...']);

% get no. of files
fileNo = length(filePaths);

% preallocate var to store epoch numbers for each file
epochNumbers = nan(fileNo, 1);

% if "loadBehav" was set to 'onesweep', preallocate a cell array for 
% storing all data in memory
if strcmp(loadBehav, 'onesweep')
    allData = cell(fileNo, 1);
end

% determine data var size in first file for the non-epoch dimensions, 
% use that in sanity checks for all other files
tmp = load(filePaths{1});
refSize = size(tmp.(varName));
% rearrange so that first dimension is epochs
epochFirstRefSize = refSize([epochDim, setdiff(1:length(refSize), epochDim)]);

% loop through files, load them
for i = 1:fileNo
    % load file content
    tmp = load(filePaths{i});
    % check size - is it what we would expect based on the first file?
    tmpSize = size(tmp.(varName));
    epochFirstTmpSize = tmpSize([epochDim, setdiff(1:length(tmpSize), epochDim)]);  % first dimension is epochs
    if ~isequal(epochFirstTmpSize(2:end), epochFirstRefSize(2:end))
        error(['Found data with unexpected size in file at ',... 
        filePaths{i}, ', investigate!']);  
    end
    % get number of epochs
    epochNumbers(i) = size(tmp.(varName), epochDim);
    % if "loadBehav" was set to 'onesweep', store data
    if strcmp(loadBehav, 'onesweep')
        allData{i} = tmp.(varName);
    end
end  % for
    
% user message
disp('Done');


%% If no epoch indices were supplied, we apply "epochMask" and / or "epochNo" if necessary.

% we only need to determine the number of epochs to include if
% "epochIndices" was not supplied
if isempty(epochIndices)

    
    %% Get indices for available epochs, apply epochMask if supplied

    % user message in case of "epochMask"
    if ~isempty(epochMask)
        disp([char(10), 'There was an epoch mask supplied, applying it ',...
            'to the epoch index lists...']);
    end

    % preallocate a cell array holding the indices of all availalbe epochs for
    % each file
    availEpochs = cell(fileNo, 1);
    % loop through files
    for i = 1:fileNo
        % if there is a non-empy "epochMask" array, we apply that to the indices
        if ~isempty(epochMask)
            % by default, all epochs are available
            tmpIndices = 1:epochNumbers(i);
            % check if the length of the mask is compatible with the number 
            % of epochs, create a temporary mask with adjusted length if
            % necessary
            if length(epochMask) < length(tmpIndices)
                tmpMask = [epochMask; zeros(length(tmpIndices)-length(epochMask), 1)];
            elseif length(epochMask) > length(tmpIndices)
                tmpMask = epochMask(1:length(tmpIndices), 1);
            else
                tmpMask = epochMask;
            end
            % apply the mask
            tmpIndices = tmpIndices(logical(tmpMask));
            % store the indices
            availEpochs{i} = tmpIndices;
            % adjust the number of available epochs
            epochNumbers(i) = numel(tmpIndices);
        % if there is no "epochMask", we just store all epochs as available
        else
            availEpochs{i} = 1:epochNumbers(i);
        end  % is ~isempty
    end  % for i

    % user message in case of "epochMask"
    if ~isempty(epochMask)
        disp('Done');
    end


    %% Determine cutoff for epoch numbers to include if "epochNo" is empty
    
    % if there was no "eopchNo" (number of epochs to include / subject)
    % supplied, get user input for it
    if isempty(epochNo)

        % report epoch numbers
        disp([char(10), 'Number of available epochs for each subject, in descending order: ']);
        disp(sort(epochNumbers, 'descend')');

        % ask for input for cutoff point
        question = [char(10), 'What should be the number of epochs kept (cutoff) for each  subject? \n',... 
            'If your choice is larger than the minimal epoch number, subjects with less epochs\n',...
            'than the cutoff will be left out of the final results array!', char(10), char(10)];
        cutoffStr = input(question, 's');

        % get numeric, check for number of remaining subjects
        epochNo = str2double(cutoffStr);
        
        % user message
        disp([char(10), 'Got it, cutoff is ', num2str(epochNo), ' epochs. ']); 
        
    end  % if isempty(epochNo)     
    
    % Check how many subjects have enough trials for the cutoff /
    % "epochNo", report to user
    if epochNo > min(epochNumbers)
        disp([char(10), 'Cutoff is larger than the minimal epoch number, ',... 
            char(10), 'data from ', num2str(sum(epochNumbers >= epochNo)),... 
            ' files out of ', num2str(fileNo), ' will remain in final data array.']);
    else
        disp([char(10), 'Cutoff is equal / smaller than the minimal epoch number, ',... 
            char(10), 'there will be epochs from all subjects in the final data array.']);
    end

   
    %% Get random epoch indices for each subject / file with enough epochs
    
    % user message
    disp([char(10), 'Generating random epoch indices for each file...']);
    
    % preallocate for selected files and corresponding epoch indices
    selectedFiles = cell(fileNo, 1);
    selectedEpochs = cell(fileNo, 1);
    
    % counter for files with enough epochs (>= epochNo)
    passCounter = 0;
    % loop through files
    for i = 1:fileNo
        % are there enough epochs?
        if epochNumbers(i) >= epochNo
            % store file path, adjust counter
            selectedFiles{i} = filePaths{i};
            passCounter = passCounter + 1;
            % get random indices
            tmpIndices = availEpochs{i};
            selectedEpochs{i} = tmpIndices(randperm(length(tmpIndices), epochNo));
        end  % if     
    end  % for i = 1:fileNo
    
    % clear empty cells
    selectedFiles(cellfun(@isempty, selectedFiles)) = [];
    selectedEpochs(cellfun(@isempty, selectedEpochs)) = [];
    
    % user message
    disp('Done');
    
    
%% If there were epoch indices supplied, perform sanity checks: 
% are they compatible with the number of files found?
% are they compatible with epoch numbers detected?

elseif ~isempty(epochIndices)
    
    % if no adjustment was needed based on "epochMask" or "epochNo", the
    % list of selected files is just the file paths cell array returned
    % earlier by listMatchingFiles
    selectedFiles = filePaths;
    
    % is the length of "epochIndices" the same as "selectedFiles"?
    if ~isequal(length(selectedFiles), length(epochIndices))
        error(['The number of epoch index lists in "epochIndices" and ',...
            'the number of detected files do not match!']);
    end
    
    % check if the supplied epoch indices are compatible with the number of
    % epochs in each file
    indexMatch = epochIdxLengths <= epochNumbers;
    if ~all(indexMatch)
        warning('More indices in "epochIndices" than available epochs for entries: ');
        disp(find(~indexMatch));
        error('The number of indices in "epochIndices" is larger than the number of available epochs for at least one subject!')
    end
    
    % if sanity checks are passed, set selectedEpochs to epochIndices
    selectedEpochs = epochIndices;
    
    % set "epochNo" from empty to the length of epoch index vectors
    epochNo = length(epochIndices{1});
    
    % user message
    disp([char(10), 'Supplied epoch indices passed minimal sanity checks, ',...
        'will use them for epoch selection.']);
    
end  % if isempty(epochIndices)


%% Loop through data, select epochs

% user message
disp([char(10), 'Selecting epochs from each file...']);

% preallocate epoch data holding var: files X epochs X [rest1 X rest2 X...]
epochData = nan([length(selectedFiles), epochNo, epochFirstRefSize(2:end)]);

% loop through data sets
for i = 1:length(selectedFiles)
    
    % if data has not been stored when epoch numbers were collected, load
    % each file in
    if strcmp(loadBehav, 'incremental')
        % get data 
        tmp = load(selectedFiles{i});
        tmpData = tmp.(varName);
        
    % else data was already stored in a cell array "allData"    
    elseif strcmp(loadBehav, 'onesweep')
        % Indexing for "allDatA" was based on the file list that might have
        % been altered since when applying "epochMask" and / or "epochNo",
        % find first the original index for current file
        tmpIdx = ismember(filePaths, selectedFiles{i});  % gives logical vector
        tmpData = allData{tmpIdx};
        
    end  % if
    
    % rearrange so that first dimension is epochs
    tmpData = permute(tmpData, [epochDim, setdiff(1:ndims(tmpData), epochDim)]);
    
    % store the epochs defined by "selectedEpochs"
    % use "hand-crafted" subscript assignment with subsref and subsasgn 
    % as we do not know for certain the dimensions of "tmpData"
    S.subs = repmat({':'}, 1, ndims(tmpData));  % subscript struct with indexing dimensions depending on "tmpData" 
    S.type = '()';  % type is a necessary field, we use the one for numeric arrays
    S.subs{1} = selectedEpochs{i};  % first dimension indices are the ones for the epochs
    tmpData = subsref(tmpData, S);  % select the relevant data slice
    % assign it to "epochData"
    S.subs = repmat({':'}, 1, ndims(epochData));
    S.subs{1} = i;
    epochData = subsasgn(epochData, S, tmpData);
           
end  % for 

% user message    
disp(['Selected epochs from all ',... 
    num2str(length(selectedFiles)), ' files.']);

% user message
disp([char(10), 'Done with everything, no errors even, this is a good day. Phew']);


return










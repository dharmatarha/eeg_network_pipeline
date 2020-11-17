function EEG = eeglabDataAggr(subFolder, subjectID, fs)
%% Function to aggregate epoch-level data (channels X samples) to 3D dataset
%
% USAGE: EEG = eeglabDataAggr(subFolder, subjectID, fs=1000)
%
% In the lab, after preprocessing in EEGLAB, we often end up with a list of
% .mat files containing epoch-level preprocessed data. This function simply
% aggregates epoch-level data for a given subject from such files into one
% 3D array (channels X samples X epochs). The returned structure is an
% EEGLAB-style struct with a few common fields defined (e.g. for subsequent
% filtering).
%
% Note that currently the function expects all epoch-level files to contain
% an "ImageGridAmp" variable holding the data for that epoch.
%
% Note the currently the function only aggregates data if epoch-level files
% all conform to the naming convention ["subjectID"_"epochNumber".mat] with
% epoch numbers forming an increasing vector of integers starting from 1.
%
% Mandatory inputs:
% subFolder - Char array, path to folder containing subject files.
% subjectID - Char array, subject identifier. Files (.mat) with epoch-level
%           data (channels X samples) are expected to start with this
%           identifier. E.g. subjectID = 'l3_s29' for files like
%           'l3_s29_100.mat'.
%
% Optional input:
% fs        - Numeric value, sampling rate in Hz. Must be integer value.
%           Defaults to 1000.
%
% Output:
% EEG       - EEGLAB-style dataset structure. EEG.data contains the
%           aggregated data (channels X samples X epochs). We set EEG.srate
%           to input arg "fs", EEG.event to empty, EEG.pnts to no. of
%           samples in one epoch, and EEG.trials to the number of epochs. 
%
% 
%


%% Input checks

if ~ismember(nargin, 2:3)
    error('Function eeglabDataAggr requires input args "subFolder" and "subjectID" while input arg "fs" is optional!');
end
if ~ischar(subFolder) || ~exist(subFolder, 'dir')
    error('Input arg "subFolder" is not a valid folder path!');
end
if ~ischar(subjectID)
    error('Input arg "subjectID" should be a character array!');
end
if nargin == 2
    fs = 1000;
else
    if ~isnumeric(fs) || length(fs)~=1 || mod(fs, 1)~=0
        error('Optional input arg "fs" must be integer value!');
    end
end

% standardize subFolder format with regards to last char
if subFolder(end) == '/'
    subFolder = subFolder(1:end-1);
end

disp([newline, 'Called function eeglabDataAggr with input args: ',...
    newline, 'folder for subject files: ', subFolder,...
    newline, 'subject identifier: ', subjectID,...
    newline, 'sampling rate: ', num2str(fs)]);


%% List all files matching subjectID in subFolder, check for epoch numbers

fileList = dir([subFolder, '/', subjectID, '*.mat']);
fileN = length(fileList);

% user feedback
disp([newline, 'Found ', num2str(fileN), ' files for subjectID ',... 
    subjectID, '. Checking if file names contain epoch numbers...']); 

% we mostly expect filenames to contain an epoch number tag after
% subjectID (e.g. "l3_s29_20.mat" where subjectID="l3_s29"), 
% check if that is the case
problemFlag = 0;
epochNumbers = nan(fileN, 1);
for i = 1:fileN
    nameParts = split(fileList(i).name(1:end-4), '_');
    if length(nameParts) == 3 && ~isempty(str2double(nameParts{end}))
        epochNumbers(i) = str2double(nameParts{end});
    else
        warning([newline, 'Cannot recognize an epoch number in the file name at ',...
            fileList(i).folder, '/', fileList(i).name,... 
            '. Cannot treat files as ordered series of epochs!']);
        problemFlag = 1;
        break;
    end
end
    
% check if epoch numbers can be sorted to a continuous ascending integer
% vector
if ~problemFlag
    if unique(diff(sort(epochNumbers))) == 1 && min(epochNumbers) == 1
        epochMin = min(epochNumbers); epochMax = max(epochNumbers);
        disp([newline, 'Epoch numbers in file names form a continuous increasing integer vector from ',... 
            num2str(epochMin), ' to ', num2str(epochMax),... 
            ', everything is as expected.']);
    else
        warning([newline, 'Epoch numbers in file names do not form a continuous ',...
            'increasing integer vector starting from 1! ',...
            'Cannot treat files as ordered series of epochs!']);
        
    end
end


%% Aggregate data from files

% only if we can sort files to epochs safely, with indices reflecting epoch
% numbers (increasing integers starting from 1)
if ~problemFlag
    
    % user feedback
    disp([newline, 'Aggregating...']);
    
    % estimate epoch data size from first file
    tmp = load([fileList(i).folder, '/', fileList(1).name]);
    [channelNo, sampleNo] = size(tmp.ImageGridAmp);
    % preallocate subject-level data var
    subData = zeros(channelNo, sampleNo, epochMax);
    
    % user feedback
    disp([newline, 'First file had data sized ', num2str(channelNo),... 
        ' x ', num2str(sampleNo), ', preallocated data var for ',... 
        num2str(epochMax), ' epochs. Loading all data...']);
    
    % loop over epochs
    for i = epochMin:epochMax
        
        % define file name for current epoch
        epochFile = [subFolder, '/', subjectID, '_', num2str(i), '.mat'];
        if ~exist(epochFile, 'file')
            error(['Could not find file for epoch no. ', num2str(i), ' at ',...
                epochFile]);
        end
        
        % load epoch-level data
        tmp = load(epochFile);
        % check its size
        if ~isequal(size(tmp.ImageGridAmp), [channelNo, sampleNo])
            error(['Unexpected data size in epoch file at ', epochFile, '!']);
        end
        
        % aggregate into 3D var
        subData(:, :, i) = tmp.ImageGridAmp;
        
    end
    
    % user feedback
    disp([newline, 'Aggregated all data']);
    
end
    

%% Generate EEGLAB-style EEG struct

if ~problemFlag
    
    % EEG struct, we define common fields used for e.g. filtering
    EEG = struct;
    EEG.data = subData;
    EEG.srate = fs;
    EEG.trials = epochMax;
    EEG.event = [];
    EEG.pnts = sampleNo;
    
    % user feedback
    disp([newline, 'Created EEGLAB-style EEG struct for data with common fields:']);
    disp(EEG);
    
end


return








function getMeanEnvelope(freq, varargin)

%% Incremental averaging of EEG envelope (or other) data.
%
% USAGE: getMeanEnvelope(freq, dirName = pwd, subjects = {'s02', 's03', ...})
%
% Calculates the group average of EEG envelope data. 
% Due to the amount of data to handle, the script uses incremental 
% averaging, calling only one dataset into memory and keeping a running 
% average.  
%
% Resulting group mean is saved into the provided
% directory (or into current working directory), named
% 'FREQUENCYBAND_envAvg.mat'.
%
% The script asssumes that the each angle .mat contains a struct "EEG" with 
% data stored in "EEG.data" (EEGlab convention, data in EEG.data). 
% Further assumption is that file naming follows the 
% 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat') convention and 
% that all relevant files are in the working directory
% 
% Mandatory input:
% freq      - Frequency band (string) to work with. Needs to be the same 
%        as used in the file names (e.g. 'alpha' for files like 
%        's01_alpha.mat'). One of {'alpha', 'beta', 'gamma', 'delta',
%        'theta'}
% 
% Optional inputs:
% dirName   - Directory path (string) pointing to the data files. Also 
%       used for saving out group means. Default is current working 
%       directory (pwd).
% subjects  - List of subjects (cell array) whose data we average. 
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

% check for mandatory argument
if nargin < 1
    error('Input arg frequency band is required!');
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end

% check optional arguments
if ~isempty(varargin)
    % if too many, raise error
    if length(varargin) > 2
        error('Too many variable inputs. Only "dirName" and "subjects" are allowed!');
    % if only one, try to find out if directory name or subjects cell array    
    elseif length(varargin) == 1
        if ischar(varargin{1})
            dirName = varargin{1};
        elseif iscell(varargin{1})
            subjects = varargin{1};
        end
        % raise error if could not identify provided argument
        if ~exist('dirName', 'var') && ~exist('subjects', 'var')
            error('At least one input arg was not recognized!');
        end     
    % if two optional arguments were provided, loop through and identify them    
    elseif length(varargin) == 2
        for i = 1:length(varargin)
            if ischar(varargin{i})
                dirName = varargin{i};
            elseif iscell(varargin{i})
                subjects = varargin{i};
            end
        end
        % raise error if at least one was not identified
        if ~exist('dirName', 'var') || ~exist('subjects', 'var')
            error('At least one input arg was not recognized!');
        end    
    end
    
end

% check if defaults are needed for input args
if ~exist('dirName', 'var')
    dirName = pwd;
end
if ~exist('subjects', 'var')
    subjects = {'s02','s03','s04','s05','s06','s07','s08','s09'...
     ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
     's21','s22','s23','s24','s25','s26','s27','s28'};
end

% user message
disp([char(10), 'Starting getMeanEnvelope function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% load first data, check dimensions
checkFile = [dirName, '/', subjects{1}, '_', freq, '.mat'];
load(checkFile, 'EEG');
[roiNo, sampleNo, epochNo, stimNo] = size(EEG.data);

% user message
disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
    num2str(stimNo), ' different stimuli (stories), ',...
    num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
    ' samples for each epoch. Assuming same dimensions for each data file.']);


%% Loop to calculate means

% user message
disp([char(10), 'Starting calculation of group mean...', char(10)])

% for creating the running average, we simply start with the first
% subject's data loaded in the previous block
runningAvg = EEG.data;

% user message
disp(['Finished with subject ', subjects{1}]);

for subIndex = 2:subNo
    
    % load next data
    nextFile = [dirName, '/', subjects{subIndex}, '_', freq,'.mat'];
    load(nextFile);    
    % check data size
    if ~isequal(size(EEG.data), [roiNo, sampleNo, epochNo, stimNo])
        error(['Data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
    end    
    % incremental averaging
    runningAvg = runningAvg + (EEG.data - runningAvg)/subIndex;

    % user message
    disp(['Finished with subject ', subjects{subIndex}]);    
    
end


%% save out the results

saveFile = [dirName, '/', freq, '_envAvg.mat'];
save(saveFile, 'runningAvg', 'subjects');

disp([char(10), 'Finished with envelope group mean', ... 
    ', saved out results to ', saveFile, char(10)]);


return

        
        
        
        
        
        
        
        
        
        
        
        
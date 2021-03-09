function [filePaths, fileNames] = listMatchingFiles(dirName, filePattern, subjects)
%% Lists files in a directory matching a pattern and a subject list
% 
% USAGE: [filePaths, fileNames] = listMatchingFiles(dirName, filePattern, subjects=[])
% 
% Util function returning both the paths and names of files from a given
% directory ("dirName") which match a given pattern ("filePatterns").
%
% Optionally it also filters for a set of subject IDs ("subjects") and
% returns only the files that containany of the subject IDs. If a subject
% list was provided, the function assumes that there is exactly one file
% per given subject, errors out otherwise.
%    
% Mandatory inputs:
% dirName           - Char array, path to folder containing files 
% filePattern       - Char array, file name part for specifying the files 
%                   to include in the selection process. The function uses
%                   "dir" for finding files, so asterisk wildcards are
%                   allowed. An extra ".mat" ending is always added when
%                   listing the files, do not include that in "filePattern"
%
% Optional input:
% subjects          - Cell array of char arrays, contains the subject IDs 
%                   for the set of subjects to include. File names are 
%                   assumed to start with the subject ID (e.g. 
%                   'l3_s10_alpha_plv.mat' for subject ID 'l3_s10'). 
%                   Determines the sequence in output arrays 
%                   "filePaths" and "fileNames". If provided, the function 
%                   assumes there is exactly one file per subject and 
%                   errors out otherwise. Defaults to [].  
%
% Outputs:
% filePaths         - Cell array of char arrays, contains the paths of
%                   files matching the inputs. If "subjects" was provided,
%                   the order of paths is determined by the IDs in
%                   subjects, so that filePaths{i} corresponds to
%                   subjects{i}.
% filenames         - Cell array of char arrays, contains the names of
%                   files matching the inputs. If "subjects" was provided,
%                   the order of paths is determined by the IDs in
%                   subjects, so that fileNames{i} corresponds to
%                   subjects{i}.
%
    

%% Input checks

if ~ismember(nargin, 2:3) 
    error(['Function listMatchingFiles requires input args "dirName", ',...
        '"filePattern" while arg "subjects" is optional!']);
end
if ~exist(dirName, 'dir')
    error('Input arg "dirName" is not a valid folder path!');
end
if ~ischar(filePattern)
    error('Input arg "filePattern" should be char array!');
end
if nargin == 2
    subjects = [];
elseif ~isempty(subjects) && (~iscell(subjects) || ~all(cellfun(@ischar, subjects)))
    error('Optional input arg "subjects" should be a cell array of char arrays or left empty!');
end


%% List matching files

% Check "dirName" and "filePattern" for data files
subFiles = dir([dirName, '/*', filePattern, '*.mat']); 
% number of files found
fileNo = length(subFiles);
% var storing file names
fileNames = extractfield(subFiles, 'name');

% user message
disp([char(10), 'Found ', num2str(fileNo),... 
    ' .mat files matching "dirName" and "filePattern"']);

% if there was a subject list supplied, check if the corresponding files
% exist, and there is only one file per cell in "subjects"
if ~isempty(subjects)
    tmpSum = zeros(size(fileNames));  % preallocate
    % go through subject ids, check if there is one file for each subject
    % id
    for s = 1:length(subjects)
        tmp = startsWith(fileNames, subjects{s});
        if sum(tmp)~=1
            error(['Error at matching data files to supplied subject ids. ',...
                'There was either multiple files or no file starting with ', subjects{s}]);
        end
        tmpSum = tmpSum+tmp;  % accumulate matches
    end
    % check if there was any file with multiple matches
    if any(tmpSum>1)
        disp(fileNames(tmpSum>1)');
        error(['Error at matching data files to supplied subject ids. ',...
            'At least one data file had more than one match to subject ids. ',...
            'See list of affected files above.']);
    end
    % if everything was fine, delete data files without a subject id match
    fileNames(tmpSum==0) = [];
    
    % force the order in "subjects" on the list of files
    % now we can take for granted that there is only one file for each
    % subject
    indices = zeros(length(subjects), 1);
    for s = 1:length(subjects)
        indices(s) = find(startsWith(fileNames, subjects{s}));
    end
    fileNames = fileNames(indices);
    
    % adjust the number of valid file names
    fileNo = length(fileNames);
    
    % user message
    disp([char(10), 'After matching potential data files to subject ids, ',...
        num2str(length(fileNames)), ' files remained.']);
end  % if

% create file paths from structs
filePaths = cell(fileNo, 1);
for i = 1:fileNo
    filePaths{i} = [dirName, '/', fileNames{i}];
end


return
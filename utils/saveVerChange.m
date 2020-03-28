function saveVerChange(filePath, newVer)
%% Change .mat version of data for python scipy compatibility
%
% USAGE: saveVersionChange(filePath, newVer='-v7')
%
% Simply loads a file a saves out its content to the same .mat file but 
% with the specified .mat version. Needed for simple load/save of .mat
% files in python with scipy.io.loadmat/savemat - they only support up to 
% version 7
%
% Mandatory input:
% filePath  - Path to .mat fle to convert.
%
% Optional input:
% newVer    - Version to convert to, defaults to "-v7"
%

%% Input checks

if ~ismember(nargin, 1:2)
    error('Function saveVerChange requires input arg "filePath" and optional arg "newVer"!');
end
if nargin==1
    newVer = '-v7';
else
    if ~ismember(newVer, {'-v4', '-v5', '-v7', '-v7.3'}) || ~ischar(newVer)
        error('Input arg "newVer" needs to be one of "-v4", "-v5", "-v7" or "-v7.3"!');
    end
end
if ~exist(filePath, 'file')
    error('Input arg "filePath" is not a valid path!');
end


%% Convert .mat file

s = load(filePath);
save(filePath, '-struct', 's', newVer);


return
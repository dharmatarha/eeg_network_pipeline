function eeglabDataAggr_wrapper(subFolder, subIds, fs)
%% Wrapper around the function eeglabDataAggr, for calling it with multiple subjects
%
% USAGE: eeglabDataAggr_wrapper(subFolder, subIds={}, fs=1000)
%
% Calls eeglabDataAggr.m for each subjectID in "subIds", passing on the
% "subFolder" and the optional "fs" args.
% The EEG structure returned by eeglabDataAggr is saved out to "subFolder"
% in a file "subjectID.mat".
%
% Mandatory inputs:
% subFolder - Char array, path to folder containing subject files.
% subIds    - Cell array holding subject identifiers. Each subject id is a
%           char array.
%           Files (.mat) with epoch-level data (channels X samples) 
%           are expected to start with these identifiers. E.g. 
%           subjectID = 'l3_s29' for files like 'l3_s29_100.mat'.
%
% Optional inputs:
% fs        - Numeric value, sampling rate in Hz. Must be integer value.
%           Defaults to 1000.
% subIds    - Cell array holding subject identifiers. Each subject id is a
%           char array.
%           Files (.mat) with epoch-level data (channels X samples) 
%           are expected to start with these identifiers. E.g. 
%           subjectID = 'l3_s29' for files like 'l3_s29_100.mat'.
%           Defaults to pre-set values for expected subject folders.


%% Input checks

% check no. of args
if ~ismember(nargin, 1:3)
    error(['Function eeglabDataAggr requires input args ',...
        '"subFolder" and "subIds" while input arg "fs" is optional!']);
end

% check mandatory arg
if ~ischar(subFolder) || ~exist(subFolder, 'dir')
    error('Input arg "subFolder" is not a valid path ot a folder!');
end

% standardize subFolder format with regards to last char
if subFolder(end) == '/'
    subFolder = subFolder(1:end-1);
end

% check optional arg fs
if nargin == 3
    if ~isnumeric(fs) || numel(fs)~=1 || mod(fs, 1)~=0
        error('Optional input arg "fs" should be integer value!');
    end
else
    fs = 1000;
end

% check optional arg subIds, set it for known folders if necessary
if ismember(nargin, 2:3)
    if ~iscell(subIds) || ~all(cellfun(@ischar, subIds))
        error('Optional input arg "subIds" should be a cell array containing char arrays!');
    end
    
else
    % assign subIds based on subject folder name 
    folderParts = split(subFolder, '/');
    lastPart = folderParts{end};
    switch lastPart
        
        case 'l3_n28'
            subIds = {'l3_s01', 'l3_s02', 'l3_s03', 'l3_s04', 'l3_s05',...
                'l3_s06', 'l3_s07', 'l3_s08', 'l3_s09', 'l3_s10',...
                'l3_s11', 'l3_s12', 'l3_s13', 'l3_s14', 'l3_s15',...
                'l3_s16', 'l3_s17', 'l3_s18', 'l3_s19', 'l3_s20',...
                'l3_s21', 'l3_s23', 'l3_s24', 'l3_s25', 'l3_s26',...
                'l3_s27', 'l3_s28', 'l3_s29'};
        case 'lh_n22'
            subIds = {'lh_s01', 'lh_s02', 'lh_s03', 'lh_s04', 'lh_s05',...
                'lh_s06', 'lh_s07', 'lh_s08', 'lh_s09', 'lh_s10',...
                'lh_s11', 'lh_s12', 'lh_s13', 'lh_s14', 'lh_s15',...
                'lh_s16', 'lh_s17', 'lh_s18', 'lh_s19', 'lh_s20',...
                'lh_s21', 'lh_s22'};
        case 'lp1_n23'
            subIds = {'lp1_s01', 'lp1_s02', 'lp1_s03', 'lp1_s04', 'lp1_s05',...
                'lp1_s07', 'lp1_s08', 'lp1_s09',...
                'lp1_s11', 'lp1_s12', 'lp1_s13', 'lp1_s14', 'lp1_s15',...
                'lp1_s16', 'lp1_s17', 'lp1_s18', 'lp1_s19', 'lp1_s20',...
                'lp1_s21', 'lp1_s22', 'lp1_s23', 'lp1_s24', 'lp1_s25'};
        case 'lp2_n25'    
            subIds = {'lp2_s01', 'lp2_s02', 'lp2_s03', 'lp2_s04', 'lp2_s05',...
                'lp2_s06', 'lp2_s07', 'lp2_s08', 'lp2_s09', 'lp2_s10',...
                'lp2_s11', 'lp2_s12', 'lp2_s13', 'lp2_s14', 'lp2_s15',...
                'lp2_s16', 'lp2_s17', 'lp2_s18', 'lp2_s19', 'lp2_s20',...
                'lp2_s21', 'lp2_s22', 'lp2_s23', 'lp2_s24', 'lp2_s25'};    
        case 'lp3_n25'    
            subIds = {'lp3_s01', 'lp3_s02', 'lp3_s03', 'lp3_s04', 'lp3_s05',...
                'lp3_s06', 'lp3_s07', 'lp3_s08', 'lp3_s09', 'lp3_s10',...
                'lp3_s11', 'lp3_s12', 'lp3_s13', 'lp3_s14', 'lp3_s15',...
                'lp3_s16', 'lp3_s17', 'lp3_s18', 'lp3_s19', 'lp3_s20',...
                'lp3_s21', 'lp3_s22', 'lp3_s23', 'lp3_s24', 'lp3_s25'};   
        case 'ls_n26'
            subIds = {'ls_s01', 'ls_s02', 'ls_s03', 'ls_s04', 'ls_s05',...
                'ls_s06', 'ls_s07', 'ls_s08', 'ls_s09', 'ls_s10',...
                'ls_s11', 'ls_s12', 'ls_s13', 'ls_s14', 'ls_s15',...
                'ls_s16', 'ls_s17', 'ls_s18', 'ls_s19', 'ls_s20',...
                'ls_s21', 'ls_s22', 'ls_s23', 'ls_s24', 'ls_s25',...
                'ls_s26'};    
        case 'lvid_n25'
            subIds = {'lvid_s01', 'lvid_s02', 'lvid_s03', 'vid_s04', 'lvid_s05',...
                'lvid_s06', 'lvid_s07', 'lvid_s08', 'lvid_s09', 'lvid_s10',...
                'lvid_s11', 'lvid_s12', 'lvid_s13', 'lvid_s14', 'lvid_s15',...
                'lvid_s16', 'lvid_s17', 'lvid_s18', 'lvid_s19', 'lvid_s20',...
                'lvid_s21', 'lvid_s22', 'lvid_s23', 'lvid_s24', 'lvid_s25'};    
        case 'lv_n28'
            subIds = {'lv_s01', 'lv_s02', 'lv_s03', 'lv_s04', 'lv_s05',...
                'lv_s06', 'lv_s07', 'lv_s08', 'lv_s09', 'lv_s10',...
                'lv_s11', 'lv_s12', 'lv_s13', 'lv_s14', 'lv_s15',...
                'lv_s16', 'lv_s17', 'lv_s18', 'lv_s19', 'lv_s20',...
                'lv_s21', 'lv_s22', 'lv_s23', 'lv_s24', 'lv_s25',...
                'lv_s26', 'lv_s27', 'lv_s28'};    
            
        otherwise
            error('Could not assign a value to optional input arg "subIds" based on arg "subFolder"!');
    end  % switch lastPart
    
end  % if nargin == 2

% user feedback
disp([newline, 'Called eeglabDataAggr_wrapper with args: ',...
    newline, 'subject folder: ', subFolder,...
    newline, 'sampling rate: ', num2str(fs),...
    newline, 'subject IDs: ']);
disp(subIds);
    

%% Loop over subjects

% number of subjects
subNo = length(subIds);

for s = 1:subNo
    
    % get current subject ID
    subject = subIds{s};
    % save file path
    saveP = [subFolder, '/', subject, '.mat'];
    
    % user feedback
    disp([newline, 'Current subject: ', subject]);
    
    % call the main function, it returns an EEGLAB-style structure
    EEG = eeglabDataAggr(subFolder, subject, fs);
    
    % save out struct
    save(saveP, 'EEG');
    
    % user feedback
    disp(['Saved out data for subject ', subject]);
    
end
    

return
    
    
    

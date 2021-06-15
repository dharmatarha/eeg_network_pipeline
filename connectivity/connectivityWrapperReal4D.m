function connectivityWrapperReal4D(freq, varargin)

%% Calculate connectivity matrices based on real valued data
%
% USAGE: connectivityWrapperReal4D(freq, dirName = pwd, subjects = {'s02', 's03', ...}, method = 'iplv')
%
% Version of connectivityWrapperReal for 4D input arrays, with dimensions
% [channels, samples, epochs, conditions].
% 
% Calculates connectivity matrices across all channels/ROIs for 
% real-valued data. Connectivity is calculated on an individual level.
%
% Connectivity is measured by calling functions defined outside this 
% script. 
% Only five connectivity measures 
% ({'plv', 'iplv', 'pli', 'ampCorr', 'orthAmpCorr'}) 
% are supported at the moment.
%
% USES PARFOR! 
%
% Results are saved into the provided directory (into current working 
% directory if omitted) separately for each subject, named
% 'SUBJECTNUMBER_FREQUENCYBAND_METHOD.mat' (e.g. 's05_alpha_iplv.mat').
%
% The function assumes that data is in EEGlab structures and that file naming
% follows the 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat')
% convention. 
% Further assumes 4D data in the structures, i.e. numel(size(EEG.data))
% = 4, with dimensions corresponding to 
% [ROI/channel, samples, epochs, conditions].
% Also assumes that all relevant files are in the working directory.
% 
% Optional input args are inferred from input arg types.
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
% subjects  - Cell array, list of subjects whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray defined in
%       "restingStateSubjects.mat" (var "subjectsRS").
% method    - String, one of {'plv', 'iplv', 'pli', 'ampCorr', 'orthAmpCorr'}. 
%       Specifies a connectivity measure. Defaults to 'iplv'.
% 


%% Input checks

% check for mandatory argument
if nargin < 1
    error('Input arg frequency band is required!');
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg freq has an unexpected value!');
end

% check optional arguments
if ~isempty(varargin) 
    disp(varargin);
    for v = 1:length(varargin)    
        if iscell(varargin{v}) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif ischar(varargin{v}) && ~exist('dirName', 'var') && exist(varargin{v}, 'dir')
            dirName = varargin{v};
        elseif ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'pli', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
            method = varargin{v};           
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "dirName", "subjects" and "method"!']);
        end
    end
end

% check if defaults are needed for input args
if ~exist('dirName', 'var')
    dirName = pwd;
end
if ~exist('subjects', 'var')
    if exist('restingStateSubjects.mat' ,'file')
        tmp = load('restingStateSubjects.mat');
        subjects = tmp.subjectsRS;
    else
        error('Couldn''t find subject list in "restingStateSubjects.mat", no subject list to work on!');
    end
end
if ~exist('method', 'var')
    method = 'iplv';
end

% user message
disp([char(10), 'Called connectivityWrapperReal4D function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% load first data, check dimensions
checkFile = [dirName, '/', subjects{1}, '_', freq, '.mat'];
load(checkFile, 'EEG');
[roiNo, sampleNo, epochNo, condNo] = size(EEG.data);
% user message
disp([char(10), 'First data file had ', num2str(condNo), ' conditions, ',... 
    num2str(epochNo), ' epochs, ', num2str(roiNo),... 
    ' channels/ROIs', ' and ', num2str(sampleNo), ... 
    ' samples. Assuming same channel & sample numbers for all subjects.']);

% prepare a lowpass filter if envelope correlation is used
if ismember(method, {'ampCorr', 'orthAmpCorr'})
    lowpassFilter = designfilt('lowpassiir',... 
            'PassbandFrequency', 10,...
            'StopbandFrequency', 13,... 
            'PassbandRipple', 0.1,... 
            'StopbandAttenuation', 50,... 
            'SampleRate', 1000,... 
            'MatchExactly',... 
            'passband');
end

% clock for timing the full function run
funcClock = tic;


%% Loop through subjects

parfor subIdx = 1:subNo
    
    % user message
    disp([char(10), 'Working on data from subject ', subjects{subIdx}]);
    
    % for timekeeping start a subject-level clock
    subClock = tic;
    
    % load subjects's data
    subFile = [dirName, '/', subjects{subIdx}, '_', freq,'.mat'];
    subData = load(subFile, 'EEG');
    subData = subData.EEG.data;
    % check data size
    [subRoiNo, subSampleNo, subEpochNo, subCondNo] = size(subData);
    if ~isequal([subRoiNo, subSampleNo], [roiNo, sampleNo])
        error(['Data for subject ', subjects{subIdx}, ' has unexpected size, investigate!']);
    end
    
    % preallocate results array
    connRes = nan(subCondNo, subEpochNo, roiNo, roiNo);
    
    % loop through conditions
    for condIdx = 1:subCondNo
    
        % loop through epochs
        for epochIdx = 1:subEpochNo

            % get phase if arg method is a phase-based measure
            if ismember(method, {'plv', 'iplv', 'pli'})
                subDataPhase = timeSeriesToPhase(squeeze(subData(:, :, epochIdx, condIdx)));  % extracts instantaneous phase from analytical signal
            end

            % measure depends on the arg method
            switch method    
                case 'plv'
                    connRes(condIdx, epochIdx, :, :) = plv(subDataPhase, 0);  % suppress messages from plv function
                case 'iplv'
                    connRes(condIdx, epochIdx, :, :) = iplv(subDataPhase, 0);  % suppress messages from iplv function    
                case 'pli'
                    connRes(condIdx, epochIdx, :, :) = pli(subDataPhase, 0);  % suppress messages from pli function
                case 'ampCorr'
                    connRes(condIdx, epochIdx, :, :) = ampCorr(squeeze(subData(:, :, epochIdx, condIdx)), lowpassFilter, 0);  % suppress messages from ampCorr function
                case 'orthAmpCorr'
                    connRes(condIdx, epochIdx, :, :) = orthAmpCorr(squeeze(subData(:, :, epochIdx, condIdx)), lowpassFilter, 0);  % suppress messages from ampCorr function            
            end  % switch method

        end  % epochIdx for loop
        
    end  % condIdx for loop
    
    % save file
    saveF = [dirName, '/' , subjects{subIdx}, '_', freq, '_', method, '.mat'];

    % save results using matfile - this is mostly allowed in a parfor loop,
    % unlike the save command
    saveM = matfile(saveF);
    saveM.connRes = connRes;
    saveM.freq = freq;
    saveM.dirName = dirName;
    saveM.method = method;
    
    % report elapsed time
    disp([char(10), 'Finished with subject ', subjects{subIdx},...
        ', it took ', num2str(toc(subClock), 3), ' secs']);
        
end  % subIdx for loop
                    

%% Cleaning up

disp([char(10), 'Done! Took ', num2str(toc(funcClock), 3),... 
    ' secs for the whole set']);


return
















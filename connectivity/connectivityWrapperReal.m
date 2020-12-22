function connectivityWrapperReal(freq, varargin)

%% Calculate connectivity matrices based on real valued data
%
% USAGE: connectivityWrapperReal(freq, dirName = pwd, subjects = {'s02', 's03', ...}, method = 'iplv')
%
% Calculates connectivity matrices across all channels/ROIs for 
% real-valued data. Connectivity is calculated on an individual level.
%
% Connectivity is measured either by calling functions (plv, iplv or pli) 
% defined outside this script or as simple amplitude correlation. 
% Only these five connectivity measures 
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
% Further assumes 3D data in the structures, i.e. numel(size(EEG.data))
% = 3, with dimensions corresponding to [ROI/channel, samples, epoch].
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
disp([char(10), 'Called connectivityWrapperReal function with following arguments: ',...
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
[roiNo, sampleNo, epochNo] = size(EEG.data);
% user message
disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
    num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
    ' samples for each epoch. Assuming same channel & sample numbers for all subjects.']);

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
    [subRoiNo, subSampleNo, subEpochNo] = size(subData);
    if ~isequal([subRoiNo, subSampleNo], [roiNo, sampleNo])
        error(['Data for subject ', subjects{subIdx}, ' has unexpected size, investigate!']);
    end
    
    % preallocate results array
    connRes = nan(subEpochNo, roiNo, roiNo);
        
    % loop through epochs
    for epochIdx = 1:subEpochNo
        
        % hacky workaround for temporary variable warnings in parfor loops
        subDataPhase = nan(roiNo);
        subDataEnv = nan(roiNo);
        
        % get phase if arg method is a phase-based measure
        if ismember(method, {'plv', 'iplv', 'pli'})
            subDataPhase = timeSeriesToPhase(squeeze(subData(:, :, epochIdx)));  % extracts instantaneous phase from analytical signal
        % get envelope for simple amplitude correlation
        elseif strcmp(method, 'ampCorr')
            subDataEnv = envelope(squeeze(subData(:, :, epochIdx))')';  % we keep the dim order channels/ROIs X samples
            % filter envelopes
            subDataEnv = filter(lowpassFilter, subDataEnv')';  % we keep the dim order channels/ROIs X samples
        end
        
        % measure depends on the arg method
        switch method    
            case 'plv'
                connRes(epochIdx, :, :) = plv(subDataPhase, 0);  % suppress messages from plv function
            case 'iplv'
                connRes(epochIdx, :, :) = iplv(subDataPhase, 0);  % suppress messages from iplv function    
            case 'pli'
                connRes(epochIdx, :, :) = pli(subDataPhase, 0);  % suppress messages from pli function
            case 'ampCorr'
                rho = corr(subDataEnv')';  % we keep the dim order channels/ROIs X samples
                rho(tril(true(roiNo))) = NaN;  % lower triangle is set to NaN, as with symmetric phase-based measures    
                connRes(epochIdx, :, :) = rho;
            case 'orthAmpCorr'
                
                % For orthogonalized amplitude correlation we use pairwise
                % orthogonalization of input time series. That is, (1) we loop
                % through all time series pairings; (2) Orthogonalize them;
                % (3) Get envelope; (4) Correlate the envelopes.
                % At the end, the upper and lower triangles of the
                % connectivity matrix are averaged out.
                
                % temp var for connectivity res
                tmp = nan(roiNo);
                
                % loops trough ROIs
                for roiOne = 1:roiNo 
                    for roiTwo = 1:roiNo
                        if roiOne ~= roiTwo
                            
                            % orthogonalize dataOne with respect to dataTwo
                            dataOne = squeeze(subData(roiOne, :, epochIdx));  
                            dataTwo = squeeze(subData(roiTwo, :, epochIdx));
                            projScalar = dot(dataOne, dataTwo)/norm(dataTwo);  % scalar projection of dataOne on dataTwo
                            projVector = projScalar*(dataTwo./norm(dataTwo));  % projection vector of dataOne on dataTwo
                            orthVector = dataOne - projVector;  % orthogonalized dataOne
                            % get envelope of both data 
                            dataOneEnv = envelope(orthVector);
                            dataTwoEnv = envelope(dataTwo);
                            % filter envelopes
                            dataOneEnv = filter(lowpassFilter, dataOneEnv);
                            dataTwoEnv = filter(lowpassFilter, dataTwoEnv);
                            % get amplitude correlation
                            tmp(roiOne, roiTwo) = corr(dataOneEnv', dataTwoEnv');
                            
                        end  % if
                    end  % for roiTwo
                end  % for roiOne
                
                % average the upper and lower triangles into the upper
                % triangle
                tmp = (triu(tmp, 1) + tril(tmp, -1)')/2;
                % lower triangle is set to NaN
                tmp(tril(true(roiNo))) = NaN;  
                % collect results in the connectivity results var
                connRes(epochIdx, :, :) = (triu(tmp, 1) + tril(tmp, -1)')/2;
                            
        end  % switch method

    end  % epochIdx for loop
    
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
















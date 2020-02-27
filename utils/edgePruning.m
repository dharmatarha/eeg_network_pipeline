function edgePruning(freq, varargin)

%% Select significant edges compared to a surrogate null model
%
% USAGE: edgePruning(freq, dirName = pwd, subjects = {'s02', 's03', ...}, method = 'plv')
%
% Calculates p-values of edges in the connectivity matrix based on
% a phase scrambling based surrogate null model.
%
% Connectivity is measured by calling functions (plv or pli) defined 
% outside this script.
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'FREQUENCYBAND_angleConnectivity.mat'.
%
% Assumes that the data is in EEGlab structures and that file naming
% follows the 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat')
% convention
% Also assumes that all relevant files are in the working directory
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
% subjects  - List of subjects (cell array) whose data we average. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
% method    - Connectivity measure compatible with phase data, one of 
%           {'plv', 'pli'}. Defaults to 'plv'.
% 
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
        elseif ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'pli', 'plv'})
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
    subjects = {'s02','s03','s04','s05','s06','s07','s08','s09'...
     ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
     's21','s22','s23','s24','s25','s26','s27','s28'};
end
if ~exist('method', 'var')
    method = 'plv';
end

% user message
disp([char(10), 'Starting angleConnectivityWrapper function with following arguments: ',...
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
[roiNo, sampleNo, epochNo, stimNo] = size(EEG.data);
% sanity check
if any(any(any(any(EEG.data>pi, 1), 2), 3), 4) || any(any(any(any(EEG.data<-pi, 1), 2), 3), 4)
    error('There is data here outside the [-pi +pi] range, is it really phase data?');
end

% user message
disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
    num2str(stimNo), ' different stimuli (stories), ',...
    num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
    ' samples for each epoch. Assuming same dimensions for each data file.']);

% Preallocate result matrix
connectivityRes = nan(subNo, stimNo, epochNo, roiNo, roiNo);
measuredConnectivityRes = nan(subNo, stimNo, epochNo, roiNo, roiNo);
prunedMeasuredConnectivityRes = nan(subNo, stimNo, epochNo, roiNo, roiNo);
meanSurrogateConnectivityRes = nan(subNo, stimNo, epochNo, roiNo, roiNo);
surrogatePValues = nan(subNo, stimNo, epochNo, roiNo, roiNo);


% save file
saveF = [dirName, '/' , freq, '_edgePruningInfo.mat'];

% clock for timing the full function run
funcClock = tic;


%% Loop through subjects

for subIdx = 1:subNo
    
    % user message
    disp([char(10), 'Working on data from subject ', subjects{subIdx}]);
    
    % for timekeeping start a subject-level clock
    subClock = tic;
    
    % load angle data - for the first subject we have already loaded it in
    % previous block
    if subIdx ~= 1%
        subFile = [dirName, '/', subjects{subIdx}, '_', freq,'.mat'];
        load(subFile, 'EEG');
        % check data size
        if ~isequal(size(EEG.data), [roiNo, sampleNo, epochNo, stimNo])
            error(['Data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
        end
        % sanity check - is it really phase data, i.e. between -pi +pi
        if any(any(any(any(EEG.data>pi, 1), 2), 3), 4) || any(any(any(any(EEG.data<-pi, 1), 2), 3), 4)
            error('There is data here outside the [-pi +pi] range, is it really phase data?');
        end
    end
    
    % loop through stimulus type / stories
    for stimIdx = 1:stimNo       
        % loop through epochs
        for epochIdx = 1:epochNo
            
            % measure depends on the arg method
            switch method    
                case 'plv'
                    eegData = squeeze(EEG.data(:, :, epochIdx, stimIdx));
                    phaseData = timeSeriesToPhase(eegData);
                    measuredConnectivityRes(subIdx, stimIdx, epochIdx, :, :) = plv(phaseData, 0);  % suppress messages from plv function
                    tic
                    surrogateConnectivityRes = generateSurrogateConnectivityData(eegData);
                    toc
                    meanSurrogateConnectivityRes(subIdx, stimIdx, epochIdx, :, :) = mean(surrogateConnectivityRes, 1);
                    surrogatePValues(subIdx, stimIdx, epochIdx, :, :) = surrogateConnectivityStatistics(squeeze(measuredConnectivityRes(subIdx, stimIdx, epochIdx, :, :)), surrogateConnectivityRes);
                    prunedMeasuredConnectivityRes(subIdx, stimIdx, epochIdx, :, :)  = pruningFunction(squeeze(surrogatePValues(subIdx, stimIdx, epochIdx, :, :)), squeeze(measuredConnectivityRes(subIdx, stimIdx, epochIdx, :, :)));

                case 'pli'
                    connectivityRes(subIdx, stimIdx, epochIdx, :, :) = pli(squeeze(EEG.data(:, :, epochIdx, stimIdx)), 0);  % suppress messages from pli function
            end
            
        end  % epochIdx for loop
    end  % stimIdx for loop
    
    % interim save
    % save(saveF, 'measuredConnectivityRes', 'meanSurrogateConnectivityRes', 'surrogatePValues', 'freq', 'dirName', 'method', 'subjects');
    
    % report elapsed time
    disp([char(10), 'Finished with subject ', subjects{subIdx},...
        ', it took ', num2str(toc(subClock), 3), ' secs']);
        
end  % subIdx for loop
                    

%% Saving, cleaning up

save(saveF, 'measuredConnectivityRes', 'prunedMeasuredConnectivityRes', 'meanSurrogateConnectivityRes', 'surrogatePValues', 'freq', 'dirName', 'method', 'subjects');

disp([char(10), 'Done! Took ', num2str(toc(funcClock), 3),... 
    ' secs for the whole set']);

return
















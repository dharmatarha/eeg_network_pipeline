function edgePruning(freq, varargin)

%% Select significant edges compared to a surrogate null model
%
% USAGE: edgePruning(freq, dirName = pwd, subjects = {'s02', 's03', ...}, method = 'plv', dRate = 1, surrNo = 10000)
%
% Calculates p-values of edges in the connectivity matrix based on
% a phase scrambling based surrogate null model.
%
% Connectivity is measured by calling functions (plv or pli) defined 
% outside this script. IMPORTANT: Only supports 'plv' as of now.
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'FREQUENCYBAND_edgePruningResults.mat'.
%
% Assumes that the data is in EEGlab structures and that file naming
% follows the 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat')
% convention
% Also assumes that all relevant files are in the working directory
% 
% Optional input args are inferred from input arg types and values.
%
% Mandatory input:
% freq      - Frequency band (string) to work with. Needs to be the same 
%        as used in the file names (e.g. 'alpha' for files like 
%        's01_alpha.mat'). One of {'alpha', 'beta', 'gamma', 'delta',
%        'theta'}
% 
% Optional inputs:
% dirName   - Directory path (string) pointing to 'angle' and 'envelope' 
%       folders. Also used for creating a 'connectivity' folder for saving 
%       out results. Default is current working directory (pwd).
% subjects  - List of subjects (cell array) whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
% method    - Connectivity measure compatible with phase data, one of 
%       {'plv', 'pli'}. Defaults to 'plv'.
% dRate     - Decimation rate. We usually work with bandpass-filtered data
%       that nevertheless retains the original (1000 or 500 Hz) sampling
%       rate. For efficiency, we can decimate this data by providing a
%       dRate integer. E.g. dRate = 10 corresponds to one 10th of the data
%       remaining after decimation. Integer between 1-20, defaults to 1 
%       (no decimation). NOTE THAT THERE IS NO ADDITIONAL FILTERING STEP 
%       INCLUDED, BEWARE OF ALIASING EFFECTS
% surrNo    - Number of surrogate data sets generated for statistical
%       testing of edge values. Integer, one of [100:100:20000], defaults 
%       to 10000. 
% 
%  


%% Input checks

% check for mandatory argument
if nargin < 1
    error('Input arg "freq" (frequency band) is required!');
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
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
        elseif isnumeric(varargin{v}) && ~exist('dRate', 'var') && ismember(varargin{v}, 1:20)
            dRate = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 100:100:20000)
            surrNo = varargin{v};
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "dirName", "subjects", "method" and "dRate"!']);
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
if ~exist('dRate', 'var')
    dRate = 1;
end
if ~exist('surrNo', 'var')
    surrNo = 10000;
end

% user message
disp([char(10), 'Starting edgePruning function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% angle and envelope folders
angleDir = [dirName, '/angle/', freq, '/'];
envDir = [dirName, '/envelope/', freq, '/'];

% add all angle and envelope files we will use to parallel loop
angleFiles = cell(subNo, 1);
envFiles = cell(subNo, 1);
for s = 1:subNo
    angleFiles{s} = [angleDir, subjects{s}, '_', freq, '.mat'];
    envFiles{s} = [envDir, subjects{s}, '_', freq, '.mat'];
end

% load first data files, both angle and envelope, check dimensions
disp([char(10), 'Sanity checks on the first angle and envelope data sets...']);
phaseData = load(angleFiles{1}, 'EEG');
phaseData = phaseData.EEG.data;
[roiNo, sampleNo, epochNo, stimNo] = size(phaseData);
% sanity check - is it phase data?
if any(any(any(any(phaseData>pi, 1), 2), 3), 4) || any(any(any(any(phaseData<-pi, 1), 2), 3), 4)
    error('There is data here outside the [-pi +pi] range, is it really phase data?');
end
% error if angle data has different size than envelope data
envData = load(envFiles{1}, 'EEG');
envData = envData.EEG.data;
if ~isequal([roiNo, sampleNo, epochNo, stimNo], size(envData))
    error('Angle and envelope data matrices have different dimensions, investigate!');
end

% decimation 
if dRate ~= 1
    phaseData = phaseData(:, 1:dRate:end, :, :);
    sampleNoOrig = sampleNo;
    sampleNo = size(phaseData, 2);
end

% user message
if dRate == 1
    disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
        ' samples for each epoch. Assuming same dimensions for each data file.']);
else
    disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNoOrig), ... 
        ' samples for each epoch. Assuming same dimensions for each data file.']);
    disp([char(10), 'After decimation (dRate = ', num2str(dRate), '), there ',...
        'are ', num2str(sampleNo), ' samples for each epoch.']);
end

% clock for timing the full function run
funcClock = tic;


%% Loop through subjects

parfor subIdx = 1:subNo
    
    % user message
    disp([char(10), 'Working on data from subject ', subjects{subIdx}]);
    
    % for timekeeping start a subject-level clock
    subClock = tic;
    
    % load angle and envelope data - for the first subject we have already 
    % loaded it in previous block, but due to how parfor works we need to
    % do so here again:
    
    % angle data
    phaseData = load(angleFiles{subIdx}, 'EEG');
    phaseData = phaseData.EEG.data;
    if dRate ~= 1
        phaseData = phaseData(:, 1:dRate:end, :, :);
    end
    % check data size
    if ~isequal(size(phaseData), [roiNo, sampleNo, epochNo, stimNo])
        error(['Angle data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
    end
    % sanity check - is it really phase data, i.e. between -pi +pi
    if any(any(any(any(phaseData>pi, 1), 2), 3), 4) || any(any(any(any(phaseData<-pi, 1), 2), 3), 4)
        error('There is data here outside the [-pi +pi] range, is it really phase data?');
    end
    % envelope data
    envData = load(envFiles{subIdx}, 'EEG');
    envData = envData.EEG.data;    
    if dRate ~= 1
        envData = envData(:, 1:dRate:end, :, :);
    end
    % check data size
    if ~isequal(size(envData), [roiNo, sampleNo, epochNo, stimNo])
        error(['Envelope data for subject ', subjects{subIdx}, ' has unexpected size, investigate!']);
    end        
    
    % preallocate result matrices
    realConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the real connectivity matrix for each epoch
    meanSurrConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the mean connectivity matrix from surrogate datasets
    pValues = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the p-values from the comparisons between real and surrogate connectivity values
    prunedConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the pruned connectivity matrices
    
    % loop through stimulus type / stories
    for stimIdx = 1:stimNo       
        % loop through epochs
        for epochIdx = 1:epochNo
            
            % get real-valued time series data from polar analytic
            [realData, ~] = pol2cart(squeeze(phaseData(:, :, epochIdx, stimIdx)), squeeze(envData(:, :, epochIdx, stimIdx))); 
            
            % measure depends on the arg method
            switch method    
                case 'plv'
                    
                    % calculate real connectivity results
                    realConn(:, :, epochIdx, stimIdx) = plv(squeeze(phaseData(:, :, epochIdx, stimIdx)), 0);
                    
                    % generate surrogate datasets and calculate
                    % connectivity matrices for them
                    surrConnData = getSurrConn(realData, surrNo);
                    % store the mean
                    meanSurrConn(:, :, epochIdx, stimIdx) = squeeze(mean(surrConnData, 1));
                    
                    % compare the real and surrogate connectivity matrices,
                    % estimate p-values for each value
                    pValues(:, :, epochIdx, stimIdx) = surrConnStats(squeeze(realConn(:, :, epochIdx, stimIdx)), surrConnData);
                    
                    % edge pruning based on pValues
                    prunedConn(:, :, epochIdx, stimIdx) = pruningFunction(squeeze(pValues(:, :, epochIdx, stimIdx)), squeeze(realConn(:, :, epochIdx, stimIdx)));

                case 'pli'
                    % NOT DONE YET
                    disp('PLI is not really supported here yet....');
%                     connectivityRes(subIdx, stimIdx, epochIdx, :, :) = pli(squeeze(EEG.data(:, :, epochIdx, stimIdx)), 0);  % suppress messages from pli function
                    
            end
            
        end  % epochIdx for loop
    end  % stimIdx for loop
    
    % get a subject-specific file name for saving results
    saveF = [dirName, '/' , subjects{subIdx}, '_', freq, '_edgePruningInfo.mat'];
    
    % save results using matfile - this is mostly allowed in a parfor loop,
    % unlike the save command
    saveM = matfile(saveF);
    saveM.realConn = realConn;
    saveM.meanSurrConn = meanSurrConn;
    saveM.pValues = pValues;
    saveM.prunedConn = prunedConn;
    saveM.method = method;
    saveM.surrNo = surrNo;
    saveM.dRate = dRate;
    saveM.subject = subjects{subIdx};
    
    % report elapsed time
    disp([char(10), 'Finished with subject ', subjects{subIdx},...
        ', it took ', num2str(toc(subClock), 3), ' secs']);
        
end  % subIdx for loop
                    

%% Saving, cleaning up

disp([char(10), 'Done! Took ', num2str(toc(funcClock), 3),... 
    ' secs for the whole set']);

return
















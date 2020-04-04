function edgePruningSingle(filePath, varargin)

%% Select significant edges compared to a surrogate null model, single dataset
%
% USAGE: edgePruningSingle(filePath, method = 'plv', dRate = 1, surrNo = 10000)
%
% Calculates p-values of edges in the connectivity matrix using
% a phase scrambling based surrogate null model.
%
% Connectivity is measured by calling functions (plv or pli) defined 
% outside this script. IMPORTANT: Only supports 'plv' as of now.
%
% Results are saved into the current working directory, named
% 'INPUT_FILENAME_edgePruningInfo.mat'.
%
% Assumes that the data is stored in a "runningAng" variable in the file.
% Data is real time series data, NOT ANGLE OR ENVELOPE DATA as expected by
% edgePruning.m
% 
% Optional input args are inferred from input arg types and values.
%
% Mandatory input:
% filePath      - String, a valid path to a .mat file. The file is expected
%           to contain a 4D "runningAvg" variable, with dimension
%           corresponding to ROIs/channels, samples, epochs and
%           stimuli/conditions
% 
% Optional inputs:
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


%% Input checks

% check for mandatory argument
if nargin < 1
    error('Input arg "filePath" is required!');
end

% check optional arguments
if ~isempty(varargin) 
    disp(varargin);
    for v = 1:length(varargin)    
        if ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'pli', 'plv'})
            method = varargin{v};     
        elseif isnumeric(varargin{v}) && ~exist('dRate', 'var') && ismember(varargin{v}, 1:20)
            dRate = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 100:100:20000)
            surrNo = varargin{v};
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "method", "dRate" and "surrNo"!']);
        end
    end
end

% check if defaults are needed for input args
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
disp([char(10), 'Starting edgePruningSingle function with following arguments: ',...
    char(10), 'File to work with: ', filePath,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo)]);


%% Basics

% load data file, check dimensions
disp([char(10), 'Loading data...']);
load(filePath, 'runningAvg');
data = runningAvg;
[roiNo, sampleNo, epochNo, stimNo] = size(data);

% decimation 
if dRate ~= 1
    data = data(:, 1:dRate:end, :, :);
    sampleNoOrig = sampleNo;
    sampleNo = size(data, 2);
end

% user message
if dRate == 1
    disp([char(10), 'Supplied data has ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
        ' samples for each epoch.']);
else
    disp([char(10), 'Supplied data has ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNoOrig), ... 
        ' samples for each epoch.']);
    disp([char(10), 'After decimation (dRate = ', num2str(dRate), '), there ',...
        'are ', num2str(sampleNo), ' samples for each epoch.']);
end

% clock for timing the full function run
funcClock = tic;


%% Calculations     
    
% preallocate result matrices
realConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the real connectivity matrix for each epoch
meanSurrConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the mean connectivity matrix from surrogate datasets
pValues = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the p-values from the comparisons between real and surrogate connectivity values
prunedConn = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the pruned connectivity matrices
    
% loop through stimulus type / stories
for stimIdx = 1:stimNo       
    % loop through epochs
    for epochIdx = 1:epochNo

        % for readability, define data for current epoch separately
        epochData = squeeze(data(:, :, epochIdx, stimIdx));
        
        % get phase data from real-valued time series
        phaseEpochData = timeSeriesToPhase(epochData); 

        % measure depends on the arg method
        switch method    
            case 'plv'

                % calculate real connectivity results
                realConn(:, :, epochIdx, stimIdx) = plv(phaseEpochData, 0);

                % generate surrogate datasets and calculate
                % connectivity matrices for them
                surrConnData = getSurrConn(epochData, surrNo);
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
    

%% Saving, return

% get a subject-specific file name for saving results
[~, fileName, ~] = fileparts(filePath);
saveF = [pwd, '/', fileName, '_edgePruningInfo.mat'];
    
% save out results
save(saveF, 'realConn', 'meanSurrConn', 'pValues', 'prunedConn', 'method',... 
    'surrNo', 'dRate', 'filePath');
    
% report elapsed time
disp([char(10), 'Finished with file ', filePath,...
    ', it took ', num2str(toc(funcClock), 3), ' secs']);


return

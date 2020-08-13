function groupLevelEdgePruning(freq, varargin)

%% Surrogate-data-based edge pruning on group-level.
%
% USAGE: groupLevelEdgePruning(freq, dirName = pwd, subjects = {'s02', 's03', ...}, surrNo = 10^4)
%
% Edge pruning on the group level, based on normal distribution fits to
% individual surrogate data sets (outputs of surrEdgeEstimation.m). 
% This method complements that implemented in edgePruning.m which prunes
% edges on the individual level, contrasting real connectivity data with 
% connectivity distributions from surrogate data.   
%
% The function relies on the outputs of surrEdgeEstimation.m that fits a
% normal or a truncated-normal to individual surrogate connectivity data.
% The parameters of the normal / truncated-normal distributions are used to
% sample from those distributions for a permutation test, comparing the 
% mean of real connectivity data with the surrogate data means, separately
% for each edge. Then an FDR (q=0.05) is applied to all edges.
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'FREQUENCYBAND_groupEdgePruningInfo.mat'.
% 
% Optional input args are inferred from input arg types and values.
%
% Mandatory input:
% freq      - Char array, one of {'alpha', 'beta', 'gamma', 'delta', 'theta'}. 
%       Frequency band to work with. Needs to be the same as used in the
%       file names (e.g. 'alpha' for files like 
%       's01_alpha_surrEdgeEstimate.mat').
% 
% Optional inputs:
% dirName   - Char array. Path of folder containing the surrEdgeEstimation 
%       output files (s*_FREQENCYBAND_surrEdgeEstimate.mat files). Also 
%       used for saving out results. Default is current working directory (pwd).
% subjects  - List of subjects (cell array) whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
% surrNo    - Number of surrogate group means generated for statistical
%       testing of edge values. Num value, one of 100:100:20000, defaults 
%       to 10^4. 
%


%% Input checks

% check for mandatory argument
if ~ismembertol(nargin, 1:4)
    error(['Function groupLevelEdgePruning requires input arg "freq" ',...
        '(frequency band) while args "dirName", "subjects" and "surrNo"',... 
        'are optional!']);
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
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 100:100:20000)
            surrNo = varargin{v};
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "dirName", "subjects" and "surrNo"!']);
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
if ~exist('surrNo', 'var')
    surrNo = 10^4;
end

% user message
disp([char(10), 'Starting groupLevelEdgePruning function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% list all surrEdgeEstimation output files we will use
for s = 1:subNo
    subjectFiles{s} = [dirName, '/', subjects{s}, '_', freq, '_surrEdgeEstimate.mat'];
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


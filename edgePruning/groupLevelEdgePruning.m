function groupLevelEdgePruning(freq, varargin)

%% Surrogate-data-based edge pruning on group-level.
%
% USAGE: groupLevelEdgePruning(realConn, freq, dirName = pwd, subjects = {'s02', 's03', ...}, surrNo = 10^4)
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
% Mandatory inputs:
% realConn  - Numeric array sized (node no. X node no. X layer no. X stim
%       no.). Contains group mean edge weights for each edge in each layer,
%       each stimulus.
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
% NOTES: 
%
%


%% Input checks

% check for mandatory argument
if ~ismembertol(nargin, 2:5)
    error(['Function groupLevelEdgePruning requires input args "realConn" ',...
        'and "freq" while args "dirName", "subjects" and "surrNo"',... 
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


%% Load and aggregate subject-level data

% number of subjects
subNo = length(subjects);

% load subject files, aggregate relevant data
for subIdx = 1:subNo
    
    % load data
    subjectFile = [dirName, '/', subjects{subIdx}, '_', freq, '_surrEdgeEstimate.mat'];
    res = load(subjectFile);
    
    % if first subject, determine dimensions and preallocate arrays holding
    % all data
    if subIdx == 1
        % check size equality of first two dimensions
        if size(res.surrNormalMu, 1)~=size(res.surrNormalMu, 2)
            error(['First two dimensions of var "surrNormalMu" in file ',... 
                subjectFile, ' should be equal (node no., node no.)!']);
        end
        % query size
        [nodeNo, ~, layerNo, stimNo] = size(res.surrNormalMu);
        
        % preallocate
        % quick check on memory requirement
        memReq = prod([nodeNo, nodeNo, layerNo, stimNo, subNo])*8*3;  % very rough estimate
        if memReq > 4*10^9
            error(['Function could use more then 4 GB memory while aggregating ',...
                'individiual data, shutting down to just to be safe...']);
        end
        % variables aggregating individual data 
        paramMu = zeros(nodeNo, nodeNo, layerNo, stimNo, subNo);
        paramSigma = paramMu;
        paramH = logical(paramMu);
        paramP = single(paramMu);  % single precision for probability values
        
        % create reference for "method" and "dRate" vars
        methodFirst = res.method;
        dRateFirst = res.dRate;
    end
    
    % check dimensions and equality across major variables
    if ~isequal(size(res.surrNormalMu), [nodeNo, nodeNo, layerNo, stimNo]) ||...
            ~isequal(size(res.surrNormalMu), size(res.surrNormalSigma)) ||...
            ~isequal(size(res.surrNormalMu), size(res.surrNormalP)) ||...
            ~isequal(size(res.surrNormalMu), size(res.surrNormalH))
        error(['At least one variable in ', subjectFile, ' has unexpected size!']);
    end
    % check "method" and "dRate" - they should be consistent across all
    % files
    if ~strcmp(res.method, methodFirst) || ~isequal(res.dRate, dRateFirst)
        error(['Variable "method" or "dRate" in ', subjectFile, ' is different than in the first file!']);
    end
    
    % aggregate data
    paramMu(:, :, :, :, subIdx) = res.surrNormalMu;
    paramSigma(:, :, :, :, subIdx) = res.surrNormalSigma;
    testH(:, :, :, :, subIdx) = logical(res.surrNormalH);
    testP(:, :, :, :, subIdx) = single(res.surrNormalP);
        
end


%% 












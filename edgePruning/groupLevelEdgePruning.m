function groupLevelEdgePruning(realConn, freq, varargin)

%% Surrogate-data-based edge pruning on group-level.
%
% USAGE: groupLevelEdgePruning(realConn, 
%                               freq, 
%                               dirName = pwd, 
%                               subjects = {'s02', 's03', ...}, 
%                               surrNo = 10^4)
%
% Edge pruning on the group level, based on (truncated) normal 
% distributions fitted to individual surrogate data sets (outputs of 
% surrEdgeEstimation.m).    
%
% The function relies on the outputs of surrEdgeEstimation.m. The function 
% surrEdgeEstimation fits a normal / truncated-normal distribution to 
% individual surrogate connectivity data.
% The parameters of normal / truncated-normal distributions are used to
% sample from those distributions for a permutation test. The test compares 
% the group-mean of real connectivity data to an ensemble of surrogate 
% data group-means. The permutation test is performed separately for each 
% edge. At the final step, an FDR (q=0.05) is applied across all edges.
%
% This method complements the one implemented in edgePruning.m which prunes
% edges on the individual level, contrasting real connectivity data with 
% connectivity distributions from surrogate data.
%
% Results are saved into the provided directory (or into current working
% directory if no "dirName" was supplied), named
% 'FREQUENCYBAND_groupEdgePruningInfo.mat'.
% 
% Optional input args are inferred from input arg types and values.
%
% Mandatory inputs:
% realConn  - Numeric array sized (node no. X node no. X layer no. X stim
%       no.). Contains group-mean edge weights (connectivity data) for 
%       each edge in each layer, each stimulus. Only upper triangles (of 
%       first two dimensions) are taken into account at calculations 
%       (assuming undirected connectivity). Might contain 0 or NaN 
%       values in upper triangles.
% freq      - Char array, one of {'alpha', 'beta', 'gamma', 'delta', 'theta'}.group stats

%       Frequency band to work with. Needs to be the same as used in the
%       file names (e.g. 'alpha' for files like 
%       's01_alpha_surrEdgeEstimate.mat').
% 
% Optional inputs:
% dirName   - Char array. Path of folder containing the surrEdgeEstimation 
%       output files (s*_FREQENCYBAND_surrEdgeEstimate.mat files). Also 
%       used for saving out results. Default is current working directory (pwd).
% subjects  - Cell array holding the list of subjects whose data we process. 
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
% Outputs: 
% Outputs are saved into a 'FREQUENCYBAND_groupEdgePruningInfo.mat' file,
% with the following variables:
% pValues   - Numeric array sized (node no. X node no. X layer no. X stim
%       no.). Contains all estimated p values from the surrogate vs. real
%       group mean comparisons, for each edge, layer/epoch  and
%       sitmulus/condition. For edges with an estimate of zero, we use the
%       value 1/surrNo*0.9, that is, the minimal possible estimate X 0.9.
% pCrit     - Numeric value. Critical p value from the FDR procedure (q =
%       .05). P values below "pCrit" are significant according to FDR. 
% fdrMask   - Boolean array sized (node no. X node no. X layer no. X stim
%       no.). Logical mask of "pValues", with true/1 where the edge
%       survived the FDR procedure ( <= pCrit).
% subjects  - Cell array, input arg "subjects".
% surrNo    - Numeric value, input arg "surrNo".
% freq      - Char array, input arg "freq".
%
%
% NOTES: 
% (1) FDR is via the utils/fdr.m function, using default options
% (Benjamini-Hochberg with q = .05). This is hardcoded behavior.
% (2) Memory requirements might become relatively large with large surrNo
% and many nodes. There is a hardcoded upper limit of ~7 GB, erroring out
% when the limit is reached. If this is a problem, either just delete the
% checks or perhaps set the var "layerSurrMeans" to single precision.
%


%% Input checks

% check no. of args
if ~ismembertol(nargin, 2:5)
    error(['Function groupLevelEdgePruning requires input args "realConn" ',...
        'and "freq" while args "dirName", "subjects" and "surrNo"',... 
        'are optional!']);
end
% check mandatory args
if ~isnumeric(realConn) || numel(size(realConn))~=4
    error('Input arg "realConn" should be a 4D numeric array!');
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
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismembertol(varargin{v}, 100:100:200000)
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
    char(10), 'Group-mean connectivity data: numeric array sized ', num2str(size(realConn)),...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'No. of surrogate group means to use for permutation tests: ', num2str(surrNo), ... 
    char(10), 'Subjects: ']);
disp(subjects);


%% Load and aggregate subject-level data

% number of subjects
subNo = length(subjects);

% bounds for truncating surrogate sample distributions
truncateBounds = [0 1];  % the phase-based connectivity values we use have range [0 1]

% user message
disp([char(10), 'Loading / aggregating all individual data...']);

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
        
        % checks on memory requirements
        % (1) rough estimate for aggregating all individual data
        dataMemReq = prod([nodeNo, nodeNo, layerNo, stimNo, subNo])*8*3; 
        % (2) rough estimate for peak memory used in inner "edgeIdx" for
        % loop (see below)
        loopMemReq = prod([surrNo, nodeNo*(nodeNo-1)/2+subNo])*8;
        % error or user message depending on memory requirements - upper bound
        % is 7 GB, so that we are ~safeish on an 8 GB machine
        disp([char(10), 'Function will attempt to use ~ ',... 
            num2str((dataMemReq+loopMemReq)/(10^9)), ' GB memory. ',...
            'This is just a rough estimate.']);
        if (dataMemReq + loopMemReq) > 7*10^9
            error(['Function could use more then 7 GB memory while aggregating ',...
                'and calculating with individiual data, shutting down to just to be safe. ',...
                'Change this hardcoded behavior or consider using single precision ',...
                'for more variables or less surrogate group means.']);
        end
        
        % preallocate
        % variables aggregating individual data 
        paramMu = zeros(nodeNo, nodeNo, layerNo, stimNo, subNo);
        paramSigma = paramMu;
        paramH = single(paramMu);  % single precision for storing NaN, 0, 1 values
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
    testH(:, :, :, :, subIdx) = single(res.surrNormalH);
    testP(:, :, :, :, subIdx) = single(res.surrNormalP);
        
end  % for subIdx

% user message
disp([char(10), 'Done loading / aggregating individual data.']);


%% Preallocate results arrays

pValues = nan(nodeNo, nodeNo, layerNo, stimNo);


%% Permutation tests

% user message
disp([char(10), 'Generating surrogate group means and comparing to real data...']);

% main clock
mainClock = tic;

% loop across stimuli / conditions
for stimIdx = 1:stimNo
    
    % user message
    disp([char(10), 'Started stimulus/condition ', num2str(stimIdx)]);
    
    % loop across layers / epochs
    for layerIdx = 1:layerNo
        
        % layer / epoch clock
        layerClock = tic;
        
        % extract layer-level group-mean connectivity values
        layerConn = squeeze(realConn(:, :, layerIdx, stimIdx));
        % linearize, only upper triangle
        layerConn(tril(true(nodeNo))) = -999;  % lower triangle is set to special value that could not occur for connectivity
        layerConn = layerConn(:);
        trilEdgeMask = layerConn==-999;  % mask of lower triangle values in linear / vector format
        layerConn(trilEdgeMask) = [];  % delete where special value
        
        % extract corresponding mu and sigma params
        layerMu = paramMu(:, :, layerIdx, stimIdx, :);
        layerMu = reshape(layerMu, [nodeNo*nodeNo, subNo]);  % one vector for a layer per subject
        layerMu(trilEdgeMask, :) = [];  % delete rows where real group-mean was deleted as lower triangle values
        layerSigma = paramSigma(:, :, layerIdx, stimIdx, :);
        layerSigma = reshape(layerSigma, [nodeNo*nodeNo, subNo]);  % one vector for a layer per subject
        layerSigma(trilEdgeMask, :) = []; % delete rows where there are only NaN values         
        
        % generate surrogate group-means "surrNo" times for each remaining edge
        
        % preallocate first
        layerSurrMeans = zeros(surrNo, nodeNo*(nodeNo-1)/2);
        % loop across edges from upper triangle
        for edgeIdx = 1:nodeNo*(nodeNo-1)/2
            
            % get truncated normal distributions for each subject, generate
            % "surrNo" samples
            
            % peallocate first
            tmpSamples = zeros(surrNo, subNo);
            % loop across subjects
            for subIdx = 1:subNo
                pd = makedist('normal', 'mu', layerMu(edgeIdx, subIdx), 'sigma', layerSigma(edgeIdx, subIdx));
                pd = truncate(pd, truncateBounds(1), truncateBounds(2));
                tmpSamples(:, subIdx) = pd.random(surrNo, 1);
            end  % for subIdx
            
            % get surrogate group-means
            layerSurrMeans(:, edgeIdx) = mean(tmpSamples, 2);
            
        end  % for edgeIdx
        
        % estimate probabilities of real group-means compared to surrogate
        % group-mean distributions
        estP = estimatedP(layerConn, layerSurrMeans, 1);  % one-tailed test  

        % p values get stored in result array with same dims as realConn
        tmp = nan(nodeNo);
        tmp(triu(true(nodeNo), 1)) = estP;
        pValues(:, :, layerIdx, stimIdx) = tmp;

        % user message
        disp([char(10), 'Done with layer/epoch ', num2str(layerIdx),... 
            ' in stimulus/condition ', num2str(stimIdx), ', took ',... 
            num2str(round(toc(layerClock), 2)), ' secs']);       
        
    end  % for layerIdx
    
end  % for stimIdx

% user message
disp([char(10), 'Done with surrogate - real group mean comparisons ',...
    'for all edges, layers/epochs and stimuli/conditions, took ',...
    num2str(round(toc(mainClock), 2)), ' secs']);


%% FDR - ignores NaN values

% user message
disp([char(10), 'FDR correction... Note that zero values in "pValues" are switched with (1/surrNo)*0.9 as a relatively conservative bound.']);

% correct for 0 values in pValues
pValues(pValues==0) = 1/surrNo*0.9;

% get critical p value at q = .05
pLin = pValues(:);
pLin(isnan(pLin)) = [];
[~, pCrit] = fdr(pLin);

% mask
fdrMask = pValues<=pCrit;

% user message
disp([char(10), 'FDR correction done, critical p is at ',... 
    num2str(pCrit), '. Saving and returning...']);


%% Save results, return

saveFile = [freq, '_groupEdgePruningInfo.mat'];
save(saveFile, 'pValues', 'fdrMask', 'pCrit', 'subjects', 'surrNo', 'freq');


return



















function simRes = connSimTest_subject_thresholding(connFile, surrDataDir, freq, method, truncated, varargin)
%% Testing the similarity of connectivity data across epoch groupings on the subject level
%
% USAGE: simRes = connSimTest_subject_thresholding(connFile, surrDataDir, freq, method, truncated, permNo=1000, metric='corr') 
%
% If given connectivity (adjacency) matrices for a set of epochs, across
% multiple subjects (in input arg "connFile"), the function 
% calculates the similarity (reliability) of average connectivity 
% matrices derived from random divisions of the epochs into two groups, 
% separately for each subject. Averaged connectivity matrices are 
% thresholded before similarity calculation (input arg "surrDataDir").
%
% Workflow, performed separately for each subject:
% (1) Randomly divide epochs into two groups
% (2) Average the connectivity (adjacency) matrices of the two groups of
% epochs
% (3) Perform thresholding of the two averaged adjacency matrices based 
% on surrogate data.
% (4) Calculate the similarity of the two averaged conenctivity matrices
% (5) Repeat steps (1)-(4) "permNo" times (default is 1000) 
%
% Current version supports similarity metrics (1) correlation (input arg
% "metric" = 'corr'), (2) Inverse Eucledian distance ('eucl'), and (3)
% DeltaCon ('deltaCon', see /networkSimilarity/deltaCon.m for details).
%
% Mandatory inputs:
% connFile              - Char array, path to group-level connectivity data
%                       file (e.g. 'group_alpha_orthAmpCorr.mat'). Contains
%                       the following variables:
%                       subjects:       cell array with subject ids
%                       epochIndices:   cell array with indices of
%                                       selected epochs for each subject
%                       connData:       numeric array containing
%                                       connectivity data, sized [subjects
%                                       X epochs X rois X rois]
% surrDataDir           - Char array, folder containing the surrogate
%                       connectivity files for each subject (e.g.
%                       'l3_s01_alpha_surrEdgeEstReal_orthAmpCorr.mat'). 
%                       Each surrogate connectivity data file should hold
%                       vars "method", "surrNormalMu", "surrNormalSigma",
%                       as usual for the outputs of the function
%                       /edgePruning/surrEdgeEstimationReal.m
% freq                  - Char array, one of {'delta', 'theta', 'alpha', 
%                       'beta', 'gamma'}. Frequency band which is reflected
%                       in surrogate connectivity file names.
% method                - Char array, one of {'plv', 'iplv', 'pli', 
%                       'ampCorr', 'orthAmpCorr'}. Connectivity method,
%                       reflected in surrogate connectivity file names.
% truncated             - Char array, one of {'truncated', 'nontruncated'}.
%                       Defines if the normal distributions derived from
%                       the surrogate data should be truncated to [0 1]. If
%                       left empty, the function tries to get this
%                       information from the surrogate data files (from var
%                       "truncated"). Defaults to empty.
%
% Optional inputs:
% truncated             - Char array, one of 
%                       {'truncated', 'nontruncated', 'fromSurrFile'}.
%                       Defines if the normal distributions derived from
%                       the surrogate data should be truncated to [0 1]. If
%                       'fromSurrFile', the function tries to get this
%                       information from the surrogate data files (from var
%                       "truncated"). Defaults to 'fromSurrFile'.
% permNo                - Numeric value, one of 10:10:10^6. Number of random
%                       permutations for epoch grouping. Defaults to 1000.
% metric                - Char array, one of {'corr', 'eucl', 'deltaCon'}.
%                       Similarity metric for comparing connectivity matrices.
%                       DeltaCon relies on the similarly named function in
%                       /networkSimilarity. Defaults to 'corr'.
%
% Output:
% simRes                - 2D numeric array sized subjects X permutations.
%                       Contains connectivity matrix similarity values for 
%                       each subject and permutation.
%


%% Input checks

% check no. of inputs
if ~ismember(nargin, 4:7)
    error(['Function connSimTest_subject_thresholding requires input args '...
    '"connFile", "surrDataDir", "freq" and "method", '...
    'while input args "truncated", "permNo" and "metric" are optional!']);
end

% check each mandatory input
if ~ischar(connFile) || ~exist(connFile, 'file')
    error(['Input arg "connFile" should be a char array, a path to a file containing ',...
        'the connectivity data from the randomly selected epochs (output of "selectEpochs")!']);
end
if ~ischar(surrDataDir) || ~exist(surrDataDir, 'dir')
    error(['Input arg "surrDataDir" should be a char array, a path to a directory containing ',...
        'the surrogate connectivity data (outputs of "surrEdgeEstimationReal")!']);
end
if ~ischar(freq) || ~ismember(freq, {'alpha', 'beta', 'gamma', 'delta', 'theta'})
    error('Input arg "freq" should be a char array, one of {"alpha", "beta", "gamma", "delta", "theta"}!');
end
if ~ischar(method) || ~ismember(method, {'plv', 'iplv', 'pli', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "method" should be a char array, one of {"plv", "iplv", "pli", "ampCorr", "orthAmpCorr"}!');
end
if ~ismember(truncated, {'truncated', 'nontruncated'})
	error('Input arg "truncated" should be one of {"truncated", "nontruncated"}!');
end

% check optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ismember(varargin{v}, {'truncated', 'nontruncated', 'fromFile'}) && ~exist('truncated', 'var')
            truncated = varargin{v};
        elseif isnumeric(varargin{v}) && ismembertol(varargin{v}, 10:10:10^6) && ~exist('permNo', 'var')
            permNo = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'corr', 'eucl', 'deltaCon'}) && ~exist('metric', 'var')
            metric = varargin{v};
        else
            error('At least one input could not be mapped nicely to args "truncated", "permNo" or "metric"!');
        end
    end
end

% assign default values if necessary
if ~exist('truncated', 'var')
    truncated = 'fromFile';
end
if ~exist('permNo', 'var')
    permNo = 1000;
end
if ~exist('metric', 'var')
    metric = 'corr';
end

% additional sanity check(s)
% check if "method" is reflected in "connFile"
if ~contains(connFile, method)
    error(['The connectivity method as defined in arg "method" is not part ',...
        'of the grouop connectivity file name ("connFile") ',...
        '- this is highly unusual. Rather safe than sorry!']);
end

% user message
disp([char(10), 'Called function sortSurrConn with inputs: ',...
    char(10), 'Connectivity file: ', connFile,...
    char(10), 'Surrogate data folder: ', surrDataDir,...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Connectivity method: ', method,...
    char(10), 'Truncate normals: ', truncated,...
    char(10), 'Number of permutations: ', num2str(permNo),...
    char(10), 'Similarity metric: ', metric]);


%% Basic params, hardcoded values

fdrQ = 0.05;  % q for FDR corrections


%% Load connectivity data

% connectivity data file should contain cell arrays for subject IDs and
% epoch indices, plus a numeric array holding all connectivity data
tmp = load(connFile);

% check names of vars in loaded connectivity data file
if ~all(ismember(fieldnames(tmp), {'subjects', 'epochIndices', 'connData'}))
    error('Connectivity data file should hold vars "subejcts", "epochIndices" and "connData"!');
end

% assign fields of interest to vars    
subjects = tmp.subjects;  % subjects X 1 cell array, each cell contains a char array as subject ID
epochIndices = tmp.epochIndices;  % subjects X 1 cell array, each cell contains a vector of epoch indices
connData = tmp.connData;  % 4D array, subects X epochs X rois X rois

% sanity checks
if size(subjects, 1) ~= size(epochIndices, 1) || size(subjects, 1) ~= size(connData, 1) 
    error('Vars "subjects", "epochIndices" and "connData" in "connFile" should have matching size on first dimension! Invalid connectivity data file?');
end
if numel(size(connData)) ~= 4
    error('Var "connData" in "connFile" should be a 4D array! Invalid connectivity data file?');
end

% get no. of subjects, epochs, ROIs from connectivity data array
[subNo, epochNo, roiNo, ~] = size(connData);

    
%% Loop through subjects

% preallocate results var
simRes = nan(subNo, permNo);
    
parfor subIdx = 1:subNo
    
    % preallocate a struct for vars in the thresholding steps and results
    
    acrossEpochs = struct;
    acrossEpochs.surrNormalMuA = nan(roiNo, roiNo);
    acrossEpochs.surrNormalMuB = nan(roiNo, roiNo);
    acrossEpochs.surrNormalSigmaA = acrossEpochs.surrNormalMuA;
    acrossEpochs.surrNormalSigmaB = acrossEpochs.surrNormalMuB;
    acrossEpochs.realConnPA = acrossEpochs.surrNormalMuA;
    acrossEpochs.realConnPB = acrossEpochs.surrNormalMuB;
    acrossEpochs.maskedConnA = acrossEpochs.surrNormalMuA;
    acrossEpochs.maskedConnB = acrossEpochs.surrNormalMuB;
    acrossEpochs.criticalPA = nan(1);
    acrossEpochs.criticalPB = nan(1);
    acrossEpochs.survivalRateA = nan(1);
    acrossEpochs.survivalRateB = nan(1);

    
    %% Load surrogate data, preparations
    
    subClock = tic;
    
    % get current subject's surrogate connectivity data file
    subID = subjects{subIdx};
    subSurrFile = dir(fullfile(surrDataDir, [subID, '_', freq, '_surrEdgeEstReal_*', method, '*.mat']));
    if numel(subSurrFile) ~= 1
        error(['There are none or too many surrogate data file(s) for subject ', subID, '!']);
    end
    subSurrPath = fullfile(subSurrFile.folder, subSurrFile.name);
    tmp = load(subSurrPath);
    
    % check if it is a file containing multiple surrogates, get method
    % index, if necessary
    methodFlag = false;
    if iscell(tmp.method) && numel(tmp.method) > 1
        methodFlag = true;
        methodIdx = find(strcmp(method, tmp.method));
        if isempty(methodIdx)
            error(['Cannot find requested method in surrogate connectivity file for subject ', subID, '!']);
        end
    end
        
    % check if arg "truncated" was set to 'fromFile', query from the loaded
    % surrogate data file if necessary
    % otherwise set subject-specific truncated value to global value
    if strcmp(truncated, 'fromFile')
        % check if there is a var in the surrogate data file names
        % "truncated"
        if ismember('truncated', fieldnames(tmp))
            % check if there should be values for multiple methods in
            % tmp.truncated
            if methodFlag
                if tmp.truncated(methodIdx)  % tmp.truncated should be logical
                    subTruncated = 'truncated';
                else
                    subTruncated = 'nontruncated';
                end
            else
                if tmp.truncated
                    subTruncated = 'truncated';
                else
                    subTruncated = 'nontruncated';
                end
            end  % if methodFlag
        
        % else: if tmp.truncated does not exist, error out    
        else
            error(['Arg "truncated" was set to "fromFile", but there is no "truncated" ',...
                'var in the surrogate connectivity data file for subject ', subID, '! ',...
                char(10), 'We cannot get the value from the file']);
        end  % if ismember  
    else
        subTruncated = truncated
    end  % if strcmp                
    
    % select mu, sigma for selected epochs,
    % change the order of their dimensions to [epochs X rois X rois],
    % as they are either [rois X rois X epochs] or [methods X rois X rois X epochs]
    if methodFlag
    	tmpMu = permute(squeeze(tmp.surrNormalMu(methodIdx, :, :, :)), [3 1 2]);
    	tmpSigma = permute(squeeze(tmp.surrNormalSigma(methodIdx, :, :, :)), [3 1 2]);     
    else
    	tmpMu = permute(tmp.surrNormalMu, [3 1 2]);
    	tmpSigma = permute(tmp.surrNormalSigma, [3 1 2]);
    end
        
    
    %% Loop through permutations
    
    for permIdx = 1:permNo
        
        % Get permutation indices, both for the connectivity array and for
        % the subject-specific surrogate data arrays
        
        % permute random epoch indices
        permIndicesA = randperm(epochNo, round(epochNo/2));
        permIndicesB = setdiff(1:epochNo, permIndicesA);
    
        % get epoch indices for current subject
        subEpochs = epochIndices{subIdx};
        if isempty(subEpochs)
            error('Epoch indices not found!');
        end
        subEpochsA = subEpochs(permIndicesA);
        subEpochsB = subEpochs(permIndicesB);

        % collect subject-level (truncated) normal params to group-level var
        surrNormalMuA = tmpMu(subEpochsA, :, :);
        surrNormalSigmaA = tmpSigma(subEpochsA, :, :);
        surrNormalMuB = tmpMu(subEpochsB, :, :);
        surrNormalSigmaB = tmpSigma(subEpochsB, :, :);

        
        %% Estimate p-values for average connectivity across epochs
        
        % Logic:
        % - Get mean connectivity across epochs
        % - For each edge, compare the real connectivity value to the normal
        % distribution derived from the per-epoch surrogate normals
        % - FDR across edges, q = 0.05
        % - Thresholding

        % mean connectivity across epochs for current subject
        subConnDataA = squeeze(mean(connData(subIdx, permIndicesA, :, :), 2));
        subConnDataB = squeeze(mean(connData(subIdx, permIndicesB, :, :), 2));

        % loops across ROIs (defining edges)
        for roi1 = 1:roiNo
            for roi2 = 1:roiNo
                
                % work only in upper triangle of connectivity matrix
                % (assume symmetric values)
                if roi1 < roi2            

                    % Derive normal distribution of surrogate data for the
                    % "group of epochs"
                    % NOTE: Linear combinations of normal distributions are
                    % just linear combinations of parameters
                    muA = mean(squeeze(surrNormalMuA(:, roi1, roi2)), 'omitnan');  % mu is simply the mean of all mu values
                    sigmaA = mean(squeeze(surrNormalSigmaA(:, roi1, roi2)), 'omitnan')/sqrt(epochNo);  % sigma is the mean sigma divided by sqrt(n)
                    muB = mean(squeeze(surrNormalMuB(:, roi1, roi2)), 'omitnan');  % mu is simply the mean of all mu values
                    sigmaB = mean(squeeze(surrNormalSigmaB(:, roi1, roi2)), 'omitnan')/sqrt(epochNo);  % sigma is the mean sigma divided by sqrt(n)

                    % Case of using simple, non-truncated normal distribution
                    if strcmp(subTruncated, 'nontruncated')
                        normalPA = normcdf(subConnDataA(roi1, roi2), muA, sigmaA);
                        normalPB = normcdf(subConnDataB(roi1, roi2), muB, sigmaB);

                    % Case of using truncated normal distribution
                    elseif strcmp(subTruncated, 'truncated')
                        % get a truncated normal distribution first
                        pdA = makedist('normal', 'mu', muA, 'sigma', sigmaA);  % returns a probability distribution object
                        pdA = truncate(pdA, 0, 1);
                        pdB = makedist('normal', 'mu', muB, 'sigma', sigmaB);  % returns a probability distribution object
                        pdB = truncate(pdB, 0, 1);
                        % get p-value from the cumulative version of the
                        % distribution function
                        normalPA = pdA.cdf(subConnDataA(roi1, roi2));
                        normalPB = pdB.cdf(subConnDataB(roi1, roi2)); 
                    end

                    % convert for a 2-two-sided test, apply correction as well
                    if normalPA > 0.5
                        normalPA = 1-normalPA;
                    end
                    if normalPB > 0.5
                        normalPB = 1-normalPB;
                    end
                    acrossEpochs.realConnPA(roi1, roi2) = normalPA*2;
                    acrossEpochs.realConnPB(roi1, roi2) = normalPB*2;

                    % store normal params
                    acrossEpochs.surrNormalMuA = muA;
                    acrossEpochs.surrNormalMuB = muB;
                    acrossEpochs.surrNormalSigmaA = sigmaA;
                    acrossEpochs.surrNormalSigmaB = sigmaB;

                end  % if roi1 < roi2 
            end  % for roi1
        end  % for roi2

        
        %% FDR corrections 
        
        % FDR correction on average connectivity matrix
        realPsA = acrossEpochs.realConnPA; 
        realPsA = realPsA(:); realPsA(isnan(realPsA)) = [];
        [~, acrossEpochs.criticalPA] = fdr(realPsA, fdrQ, 'bh');
        realPsB = acrossEpochs.realConnPB; 
        realPsB = realPsB(:); realPsB(isnan(realPsB)) = [];
        [~, acrossEpochs.criticalPB] = fdr(realPsB, fdrQ, 'bh'); 

        % create masked connectivity tensor
        pMaskA = acrossEpochs.realConnPA <= acrossEpochs.criticalPA;
        tmpDataA = subConnDataA;
        tmpDataA(~pMaskA) = 0;  % values not surviving the threshold are set to zero
        acrossEpochs.maskedConnA = tmpDataA;
        pMaskB = acrossEpochs.realConnPB <= acrossEpochs.criticalPB;
        tmpDataB = subConnDataB;
        tmpDataB(~pMaskB) = 0;  % values not surviving the threshold are set to zero
        acrossEpochs.maskedConnB = tmpDataB;

        % get rate of edges surviving pruning
        acrossEpochs.survivalRateA = sum(pMaskA(:), 'omitnan')/(roiNo*(roiNo-1)/2);
        acrossEpochs.survivalRateB = sum(pMaskB(:), 'omitnan')/(roiNo*(roiNo-1)/2);
        
        
        %% calculate similarity of thresholded matrices
        
        matrixA = acrossEpochs.maskedConnA;
        matrixB = acrossEpochs.maskedConnB;
        
        % for certain metrics get symmetric adjacency matrices with zeros at diagonal
        if ismember(metric, {'eucl', 'deltaCon'})
            matrixA = triu(matrixA, 1) + triu(matrixA, 1)'; 
            matrixB = triu(matrixB, 1) + triu(matrixB, 1)';
            % also standardize the scale of connections across the two
            % matrices to a common sum (=10)
            matrixA = 10*matrixA./sum(matrixA(:));
            matrixB = 10*matrixB./sum(matrixB(:));
        end

        % calculate similarity according to arg "metric"
        switch metric
            
            case 'corr'
                % linearize upper triangles of mean connectivity matrices
                linA = matrixA(triu(true(roiNo), 1));
                linB = matrixB(triu(true(roiNo), 1));
                simRes(subIdx, permIdx) = corr(linA, linB);
                
            case 'eucl'
                % we use a distance measure if 'eucl' is selected
                % (Frobenius norm)
                simRes(subIdx, permIdx) = norm(matrixA-matrixB, 'fro');
                
            case 'deltaCon'
                % we use a distance measure (DeltaCon distance) if
                % 'deltaCon' is selected (second output of deltaCon.m)
                [~, simRes(subIdx, permIdx)] = deltaCon(matrixA, matrixB, false);  % "false" is for verbosity
                
        end  % switch metric
        
    end % for permIdx
    
    
    % display elapsed time 
    elapsedT = round(toc(subClock), 2);
    disp([char(10), 'Finished with subject ', subID, '. Took ', num2str(elapsedT), ' secs.'])

    
end % for subIdx

return
        
        













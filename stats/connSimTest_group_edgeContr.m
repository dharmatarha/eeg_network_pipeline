function [simRes, edgeContr] = connSimTest_group_edgeContr(connArray, varargin)
%% Testing the similarity of connectivity data across epoch groupings on the group level and the contribution of each edge to the similarity
%
% USAGE: simRes = connSimTest_group_edgeContr(connArray, permNo=1000, metric='corr') 
%
% If given connectivity (adjacency) matrices for a set of epochs, across
% multiple subjects (in input arg "connArray"), the function calculates the
% similarity (reliability) of average connectivity matrices derived from 
% random divisions of the dataset into two groups (keeping all epochs for
% each subject within a single group). The function also calculates the
% contribution of each edge to the similarity.
%
% Workflow, separately performed for each subject:
% (1) Randomly divide epochs into two groups (keeping all epochs for
% each subject within a single group)
% (2) Average the connectivity (adjacency) matrices of the two groups of
% epochs
% (3) Calculate the similarity of the two averaged conenctivity matrices
% and the contribution of each edge to the similarity.
% (4) Repeat steps (1)-(3) "permNo" times (default is 1000) 
%
% At the moment, support similarity metrics (1) correlation (input arg
% "metric" = 'corr'), (2) Inverse Eucledian distance ('eucl'), and (3)
% DeltaCon ('deltaCon', see /networkSimilarity/deltaCon.m for details).
%
% Mandatory inputs:
% connArray     - 4D numeric array, with dimensions: subjects X epochs X
%               ROIs X ROIs. Contains a connectivity matrix for each epoch
%               of each subject. Might only contain valid data in the upper
%               triangles of the connectivity matrices, the rest might be
%               NaN (in case of symmetric connectivity measures).
%
% Optional inputs:
% permNo        - Numeric value, one of 10:10:10^6. Number of random
%               permutations for epoch grouping. Defaults to 1000.
% metric        - Char array, one of {'corr', 'eucl', 'deltaCon'}.
%               Similarity metric for comparing connectivity matrices.
%               DeltaCon relies on the similarly named function in
%               /networkSimilarity. Defautls to 'corr'.
%
% Output:
% simRes        - 1D numeric array sized 1 X permutations. Contains
%               group similarity values for each permutation.
%
% edgeContr     - 2D numeric array sized (number of edges) X permutations.
%               Contains edge contribution values to the similarity
%               for each permutation and edge.
%


%% Input checks

% cehck no. if inputs
if ~ismember(nargin, 1:3)
    error('Function connSimTest_group requires input arg "connArray" while args "permNo" and "metric" are optional!');
end
% check mandatory input
if ~isnumeric(connArray) || numel(size(connArray))~=4 || size(connArray, 3)~=size(connArray ,4)
    error('Input rag "connArray" should be a 4D numeric array where the 3rd and 4th dimensions have the same size!');
end
% check optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismembertol(varargin{v}, 10:10:10^6) && ~exist('permNo', 'var')
            permNo = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'corr', 'eucl', 'deltaCon'}) && ~exist('metric', 'var')
            metric = varargin{v};
        else
            error('At least one input could not ba mapped nicely to args "permNo" or "metric"!');
        end
    end
end
% assign default values if necessary
if ~exist('permNo', 'var')
    permNo = 1000;
end
if ~exist('metric', 'var')
    metric = 'corr';
end

% user message
disp([char(10), 'Called connSimTest_group with input args: ',...
    char(10), 'Input array sized ', num2str(size(connArray)), ...
    char(10), 'Numer of permutations: ', num2str(permNo),...
    char(10), 'Similarity metric: ', metric]);


%% Calculate similarities

% clock for whole run
funcClock = tic;

% get subject, epoch and channel/ROI numbers
[subNo, ~, roiNo, ~] = size(connArray);

% preallocate results var
simRes = nan(1, permNo);
edgeContr = nan(((roiNo*roiNo-roiNo)/2), permNo);

% average over epochs before subject-to-subject
% comparisons
connArray = squeeze(mean(connArray, 2));

% user message
disp([char(10), 'Calculating...']);

% loop through permutations
for permIdx = 1:permNo
    
    % permute random subject indices
    permIndicesA = randperm(subNo, round(subNo/2));
    permIndicesB = setdiff(1:subNo, permIndicesA);
    
    % divide subjects according to random indices 
    permTensorA = squeeze(connArray(permIndicesA, :, :));
    permTensorB = squeeze(connArray(permIndicesB, :, :));
    
    % average subjects' data in the two groups
    permMatrixA = squeeze(mean(permTensorA, 1));
    permMatrixB = squeeze(mean(permTensorB, 1));    
    
    % for certain metrics get symmetric adjacency matrices with zeros at diagonal
    if ismember(metric, {'eucl', 'deltaCon'})
        permMatrixA = triu(permMatrixA, 1) + triu(permMatrixA, 1)';
        permMatrixB = triu(permMatrixB, 1) + triu(permMatrixB, 1)';
        % also standardize the scale of connections across the two
        % matrices to a common sum (=10)
        permMatrixA = 10*permMatrixA./sum(permMatrixA(:));
        permMatrixB = 10*permMatrixB./sum(permMatrixB(:));        
    end
    
    % calculate similarity according to arg "metric"
    switch metric
        
        case 'corr'
            % linearize upper triangles of mean connectivity matrices
            linA = permMatrixA(triu(true(roiNo), 1));
            linB = permMatrixB(triu(true(roiNo), 1));
            simRes(permIdx) = corr(linA, linB);
            
            zA = (linA - mean(linA)) ./ norm(linA - mean(linA));
            zB = (linB - mean(linB)) ./ norm(linB - mean(linB));
            edgeContr(:, permIdx) = (zA .* zB) ./ sum(zA .* zB);
            
        case 'eucl'
%             % similarity is based on norm of difference
%             simRes(permIdx) = 1/(1 + norm(permMatrixA-permMatrixB, 'fro'));
            % we use a distance measure if 'eucl' is selected
            % (Frobenius norm)
            simRes(permIdx) = norm(permMatrixA-permMatrixB, 'fro');            
            
        case 'deltaCon'
%             % first output arg of deltaCon is similarity
%             simRes(permIdx) = deltaCon(permMatrixA, permMatrixB, false);  % "false" is for verbosity
            % we use a distance measure (DeltaCon distance) if
            % 'deltaCon' is selected (second output of deltaCon.m)
            [~, simRes(permIdx)] = deltaCon(permMatrixA, permMatrixB, false);  % "false" is for verbosity            
            
    end  % switch metric
    
end  % for permIdx


%% End message, return

% user message
disp([char(10), 'Done with everything, running the function took ',... 
    num2str(round(toc(funcClock), 2)), ' secs']);


return
        
        













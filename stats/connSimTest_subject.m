function simRes = connSimTest_subject(connArray, varargin)
%% Testing the similarity of connectivity data across epoch groupings on the subject level
%
% USAGE: simRes = connSimTest_subject(connArray, permNo=1000, metric='corr') 
%
% Mandatory inputs:
% connArray     - 4D numeric array, with dimensions subjects X epochs X
%               ROIs X ROIs. Contains a connectivity matrix for each epoch
%               of each subject. Might only contain valid data in the upper
%               triangles of the connectivity matrices, the rest might be
%               NaN (in case of symmetric connectivity measures).
% Optional inputs:
% permNo        - Numeric value, one of 10:10:10^6. Number of random
%               permutations for epoch grouping. Defaults to 1000.
% metric        - Char array, one of {'corr', 'eucl', 'deltaCon'}.
%               Similarity metric for comparing connectivity matrices.
%               DeltaCon relies on the similarly named function in
%               /networkSimilarity. Defautls to 'corr'.
%
% Output:
% simRes        - 2D numeric array sized subjects X permutations. Contains
%               connectivity matrix similarity values for each subject and
%               permutation.
%


%% Input checks

if ~ismember(nargin, 1:3)
    error('Function connSimTest_subject requires input arg "connArray" while args "permNo" and "metric" are optional!');
end
if ~isnumeric(connArray) || numel(size(connArray))~=4 || size(connArray, 3)~=size(connArray ,4)
    error('Input rag "connArray" should be a 4D numeric array where the 3rd and 4th dimensions have the same size!');
end
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
if ~exist('permNo', 'var')
    permNo = 1000;
end
if ~exist('metric', 'var')
    metric = 'corr';
end

disp([char(10), 'Called connSimTest_subject with input args: ',...
    char(10), 'Input array sized ', num2str(size(connArray)), ...
    char(10), 'Numer of permutations: ', num2str(permNo),...
    char(10), 'Similarity metric: ', metric]);



%% Calculate similarities

% clock for whole run
funcClock = tic;

% get subject, epoch and channel/ROI numbers
[subNo, epochNo, roiNo, ~] = size(connArray);

% preallocate results var
simRes = nan(200, permNo);

% user message
disp([char(10), 'Calculating...']);

% loop through subjects
for subIdx = 1:subNo
    
    % subject-level clock
    subClock = tic;
    
    % select subject's data
    subData = squeeze(connArray(subIdx, :, :, :));
    
    % loop through permutations
    for permIdx = 1:permNo
        
        % permute random epoch indices
        permIndicesA = randperm(epochNo, round(epochNo/2));
        permIndicesB = setdiff(1:epochNo, permIndicesA);
        % divide epochs according to random indices and average them in two
        % groups
        permMatrixA = squeeze(mean(subData(permIndicesA, :, :), 1));
        permMatrixB = squeeze(mean(subData(permIndicesB, :, :), 1));
        
        % for certain metrics get symmetric adjacency matrices with zeros at diagonal
        if ismember(metric, {'eucl', 'deltaCon'})
            permMatrixA = triu(permMatrixA, 1) + triu(permMatrixA, 1)'; 
            permMatrixB = triu(permMatrixB, 1) + triu(permMatrixB, 1)';
        end

        % calculate similarity according to arg "metric"
        switch metric
            
            case 'corr'
                % linearize upper triangles of mean connectivity matrices
                linA = permMatrixA(triu(true(roiNo), 1));
                linB = permMatrixB(triu(true(roiNo), 1));
                simRes(subIdx, permIdx) = corr(linA, linB);
                
            case 'eucl'
                % similarity is based on norm of difference
                simRes(subIdx, permIdx) = 1/(1 + norm(permMatrixA-permMatrixB, 'fro'));
                
            case 'deltaCon'
                % first output arg of deltaCon is similarity
                simRes(subIdx, permIdx) = deltaCon(permMatrixA, permMatrixB, false);  % "false" is for verbosity
                
        end  % switch metric
        
    end  % for permIdx
    
    % user message
    disp([char(10), 'Done with subject ', num2str(subIdx),... 
        ', took ', num2str(round(toc(subClock), 2)), ' secs']);
    
end  % for subIdx


%% End message, return

% user message
disp([char(10), 'Done with everything, running the function took ',... 
    num2str(round(toc(funcClock), 2)), ' secs']);


return
        
        













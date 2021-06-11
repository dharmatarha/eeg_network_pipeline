function [distRes] = connDistanceTest_betweenSubject_epochAveraged(connArray, varargin)
%% Testing the distance of connectivity data between different subjects
%
% USAGE: distRes = connDistanceTest_betweenSubject_epochAveraged(connArray, metric='corr') 
%
% If given connectivity (adjacency) matrices (epoch-averaged) across
% multiple subjects (in input arg "connArray"), the function calculates 
% the distance of connectivity matrices between different subjects.
%
% Workflow:
% Calculate the distance of two averaged conenctivity matrices
% for all possible pairing of all subjects.
%
% At the moment, support distance metrics (1) correlation (input arg
% "metric" = 'corr'), (2) Eucledian distance ('eucl'), (3) adjacency spectral 
% distance ('adjacencySpectral'), (4) Laplacian spectral distance ('LaplacianSpectral')
% and (5) DeltaCon ('deltaCon', see /networkSimilarity/deltaCon.m for details).
%
% Mandatory inputs:
% connArray     - 3D numeric array, with dimensions: subjects X ROIs X ROIs.
%               Contains an epoch-averaged connectivity matrix for each subject.
%               Might only contain valid data in the upper triangles of the 
%               connectivity matrices, the rest might be NaN (in case of symmetric
%               connectivity measures).
%
% Optional inputs:
% metric        - Char array, one of {'corr', 'eucl', 'adjacencySpectral', 
%               'LaplacianSpectral', 'deltaCon'}.
%               Distance metric for comparing connectivity matrices.
%               adjacencySpectral, LaplacianSpectral and DeltaCon rely on 
%               similarly named functions in
%               /networkSimilarity. Defautls to 'corr'.
%
% Output:
% distRes        - 2D numeric array sized subjects X subjects. Contains
%               connectivity matrix distance values for all possible
%               pairings of all subjects. Only upper triangle is populated.
%
% Note: Correlation is basically a similarity metric. It is converted to distance
%       using the formula: distance = sqrt(1 - similarity).
%       Negative correlation values are transformed to non-negative values
%       by a scale transformation.

%% Input checks

% check no. if inputs
if ~ismember(nargin, 1:2)
    error('Function connDistanceTest_betweenSubject requires input arg "connArray" while arg "metric" is optional!');
end
% check mandatory input
if ~isnumeric(connArray) || numel(size(connArray))~=3 || size(connArray, 2)~=size(connArray ,3)
    error('Input arg "connArray" should be a 3D numeric array where the 2nd and 3rd dimensions have the same size!');
end
% check optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ismember(varargin{v}, {'corr', 'eucl', 'adjacencySpectral', 'LaplacianSpectral', 'deltaCon'}) && ~exist('metric', 'var')
            metric = varargin{v};
        else
            error('Input could not be mapped nicely to arg "metric"!');
        end
    end
end
% assign default values if necessary
if ~exist('metric', 'var')
    metric = 'corr';
end

% user message
disp([char(10), 'Called connDistanceTest_betweenSubject with input args: ',...
    char(10), 'Input array sized ', num2str(size(connArray)), ...
    char(10), 'Distance metric: ', metric]);


%% Calculate similarities

% clock for whole run
funcClock = tic;

% get subject and channel/ROI numbers
[subNo, roiNo, ~] = size(connArray);

% preallocate results var
distRes = nan(subNo, subNo);

% user message
disp([char(10), 'Calculating...']);

% loop through subjects
for subIdx = 1:subNo
    
    % subject-level clock
    subClock = tic;
    
    % select subject's connectivity data (averaged over epochs)
    subMatrixA = squeeze(connArray(subIdx, :, :));
    
    % get symmetric adjacency matrices with zeros at diagonal
    subMatrixA = triu(subMatrixA, 1) + triu(subMatrixA, 1)';
    % also standardize the scale of connections across the two
    % matrices to a common sum (=10)
    subMatrixA = 10*subMatrixA./sum(subMatrixA(:));
    
    % loop through compared subjects
    for compSubIdx = subIdx+1 : subNo
        
        % select compared subject's connectivity data (averaged over
        % epochs)
        subMatrixB = squeeze(connArray(compSubIdx, :, :));
        
        % get symmetric adjacency matrices with zeros at diagonal
        subMatrixB = triu(subMatrixB, 1) + triu(subMatrixB, 1)';
        % also standardize the scale of connections across the two
        % matrices to a common sum (=10)
        subMatrixB = 10*subMatrixB./sum(subMatrixB(:));

        % calculate similarity according to arg "metric"
        switch metric
            
            case 'corr'
                % linearize upper triangles of mean connectivity matrices
                linA = subMatrixA(triu(true(roiNo), 1));
                linB = subMatrixB(triu(true(roiNo), 1));
                distRes(subIdx, compSubIdx) = corr(linA, linB);              
                
            case 'eucl'
%                 % similarity is based on norm of difference
%                 simRes(subIdx, compSubIdx) = 1/(1 + norm(subMatrixA-subMatrixB, 'fro'));
                % we use a distance measure if 'eucl' is selected
                % (Frobenius norm)
                distRes(subIdx, compSubIdx) = norm(subMatrixA-subMatrixB, 'fro');
                
            case 'adjacencySpectral'
                % first output arg of adjacencySpectralDistance is similarity
                [~, distRes(subIdx, compSubIdx)] = adjacencySpectralDistance(subMatrixA, subMatrixB, false);
                
            case 'LaplacianSpectral'
                % first output arg of laplacianSpectralDistance is similarity
                [~, distRes(subIdx, compSubIdx)] = laplacianSpectralDistance(subMatrixA, subMatrixB, false);
                
            case 'deltaCon'
%                 % first output arg of deltaCon is similarity
%                 simRes(subIdx, compSubIdx) = deltaCon(subMatrixA, subMatrixB, false);  % "false" is for verbosity
                % we use a distance measure (DeltaCon distance) if
                % 'deltaCon' is selected (second output of deltaCon.m)
                [~, distRes(subIdx, compSubIdx)] = deltaCon(subMatrixA, subMatrixB, false);  % "false" is for verbosity                
                
        end  % switch metric
        
    end  % for compSubIdx
    
    % user message
    disp([char(10), 'Done with subject ', num2str(subIdx),... 
        ', took ', num2str(round(toc(subClock), 2)), ' secs']);
    
end  % for subIdx
        
if ismember(metric, {'corr'})
    smallestElement = min(min(distRes));
    distRes = (distRes - smallestElement) ./ (1-smallestElement);
    distRes = sqrt(1 - distRes);
end

%% End message, return

% user message
disp([char(10), 'Done with everything, running the function took ',... 
    num2str(round(toc(funcClock), 2)), ' secs']);


return
        
        













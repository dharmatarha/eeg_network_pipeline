function [similarity, distance] = adjacencySpectralDistance(adjMatrix1, adjMatrix2, varargin)
%% Calculate adjacency spectral distance between two graphs 
%
% USAGE: [similarity, distance] = adjacencySpectralDistance(adjMatrix1, adjMatrix2, k = N, verbose = true)
% 
% The function calculates the adjacency spectral distance for the
% two networks defined by the adjancency matrices "adjMatrix1" and
% "adjMatrix2".  See the details in:
% Wilson et al. (2008). A study of graph spectra for comparing graphs and trees
%   
% Briefly, the measure calculates the eigenvalues of the input matrices 
% and calculates the distance between the two spectra in the l2 metric.
%
% Does not require node correspondence, works on both
% weighted and unweighted graphs.
%
% Mandatory inputs:
% adjMatrix1    - Numeric square matrix, adjacency matrix for the first 
%               network. Values on the diagonal must be zeros.
% adjMatrix2    - Numeric square matrix, adjacency matrix for the second 
%               network, same size as "adjMatrix1". Values on the diagonal 
%               must be zeros.
%
% Optional input:
% k             - Numeric value, number of the largest eigenvalues to
%               compare. Defaults to the number of nodes.
% verbose       - Logical value (true or false). Verbosity, "false" 
%               meaning no user messages, "true" meaning user messages.
%
% Outputs:
% similarity    - Numeric value, adjacency spectral similarity.
% distance      - Numeric value, adjacency spectral distance.
%


%% Input checks

% check number of arguments
if ~ismember(nargin, 2:4)
    error(['Function adjacencySpectralDistance requires input args "adjMatrix1" and "adjMatrix2", ',...
        'while input args "k" and "verbose" are optional!']);
end
% check mandatory args
if ~isnumeric(adjMatrix1) || ~ismatrix(adjMatrix1) || size(adjMatrix1, 1)~=size(adjMatrix1, 2)
    error(['Input arg "adjMatrix1" should be a numeric square matrix ',...
        '(adjacency matrix of a network)!']);
end
if ~isnumeric(adjMatrix2) || ~ismatrix(adjMatrix2) || size(adjMatrix2, 1)~=size(adjMatrix2, 2)
    error(['Input arg "adjMatrix2" should be a numeric square matrix ',...
        '(adjacency matrix of a network)!']);
end
if ~isequal(size(adjMatrix1), size(adjMatrix2))
    error('Input args "adjMatrix1" and "adjMatrix2" should have the same size!');
end
% check optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('k', 'var')
            k = varargin{v};
        elseif islogical(varargin{v}) && numel(varargin{v})==1 && ~exist('verbose', 'var')
            verbose = varargin{v};
        else
            error('An optional input arg does not match nicely to "k" or "verbose"!');
        end
    end
end
% assign defaults
if ~exist('k', 'var')
    k = size(adjMatrix1, 1);
end
if ~exist('verbose', 'var')
    verbose = true;
end
% any further check
if any(diag(adjMatrix1)) || any(diag(adjMatrix2))
    error('There is at least one non-zero value on one of the diagonals!');
end
if k > size(adjMatrix1, 1)
    k = size(adjMatrix1, 1);
    warning('k was limited to the number of nodes');
end

% user message if verbose
if verbose
    disp([char(10), 'Called adjacencySpectralDistance function with input args: ',...
        char(10), 'Adjacency (connectivity) matrices of size ', num2str(size(adjMatrix1)),...
        char(10), 'k (number of the largest eigenvalues to compare): ', num2str(k),...
        char(10), 'Verbosity: ', num2str(verbose), char(10)]);
end


%% Calculate adjacency spectral distance

% number of nodes
nodeNo = size(adjMatrix1, 1);

% adjacency spectra
lambda1 = eig(adjMatrix1);
lambda2 = eig(adjMatrix2);

% keeping the largest k eigenvalues
lambda1 = sort(lambda1, 'descend');
lambda1 = lambda1(1:k);
lambda2 = sort(lambda2, 'descend');
lambda2 = lambda2(1:k);

% distance
distance = sqrt(sum((lambda1 - lambda2).^2));

% similarity is bounded to [0 1]
similarity = 1/(1+distance);

% user message if verbose
if verbose
    disp([char(10), 'Adjacency spectral distance: ', num2str(distance),...
        char(10), 'Adjacency spectral similarity: ', num2str(similarity), char(10)]);
end


return
















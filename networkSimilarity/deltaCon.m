function [deltaConSim, deltaConDist, S1, S2] = deltaCon(adjMatrix1, adjMatrix2, varargin)
%% Calculate DeltaCon network similarity measure 
%
% USAGE: [deltaConSim, deltaConDist, S1, S2] = deltaCon(adjMatrix1, adjMatrix2, epsilon = 0.01, verbose = true)
% 
% The function calculates the DeltaCon network similarity measure for the
% two networks defined by the adjancency matrices "adjMatrix1" and
% "adjMatrix2". We implement the slower, exact algorithm that has 
% quadratic complexity in node number. See the details in:
%   Koutra et al. (2013). Deltacon: A principled massive-graph 
%   similarity function. 
%   Koutra et al. (2016). DELTACON: Principled Massive-Graph Similarity
%   Function with Attribution
%   
% Briefly, the measure first calculates pairwise node affinities for both 
% input matrices using Fast Belief Propagation and then compares the node
% affinity matrices with root euclidean distance (a.k.a. Matusita dist). 
%
% Requires identical node sets across the two networks, works on both
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
% epsilon       - Numeric value, small constant. Controls the relative 
%               weight of node ranks and the adjacency matrix. See the 
%               paper and the equations for details. Defaults to 0.01.
% verbose       - Logical value (true or false). Verbosity, "false" 
%               meaning no user messages, "true" meaning user messages.
%
% Outputs:
% deltaConSim   - Numeric value, DeltaCon similarity.
% deltaConDist  - Numeric value, DeltaCon distance. 
% S1            - Numeric matrix, pairwise node affinities for
%               "adjMatrix1". Same size as "adjMatrix1" (nodes X nodes).
% S2            - Numeric matrix, pairwise node affinities for
%               "adjMatrix2". Same size as "adjMatrix2" (nodes X nodes).
%
% NOTES:
% (1) Adjust epsilon if deltaConSim is complex
%


%% Input checks

% check number of arguments
if ~ismember(nargin, 2:4)
    error(['Function deltaCon requires input args "adjMatrix1" and "adjMatrix2", ',...
        'while input args "epsilon" and "verbose" are optional!']);
end
% check mandatory args
if ~isnumeric(adjMatrix1) || ~ismatrix(adjMatrix1) || size(adjMatrix1, 1)~=size(adjMatrix1, 2)
    error(['Input arg "adjMatrix1" should be a numeric square matrix ',...
        '(adjacency matrix of a network)!']);
end
if ~isnumeric(adjMatrix2) || ~ismatrix(adjMatrix2) || size(adjMatrix2, 1)~=size(adjMatrix2, 2)
    error(['Input arg "adjMatrix1" should be a numeric square matrix ',...
        '(adjacency matrix of a network)!']);
end
if ~isequal(size(adjMatrix1), size(adjMatrix2))
    error('Input args "adjMatrix1" and "adjMatrix2" should have the same size!');
end
% check optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('epsilon', 'var')
            epsilon = varargin{v};
        elseif islogical(varargin{v}) && numel(varargin{v})==1 && ~exist('verbose', 'var')
            verbose = varargin{v};
        else
            error('An optional input arg does not match nicely to "epsilon" or "verbose"!');
        end
    end
end
% assign defaults
if ~exist('epsilon', 'var')
    epsilon = 0.01;
end
if ~exist('verbose', 'var')
    verbose = true;
end
% any further check
if any(diag(adjMatrix1)) || any(diag(adjMatrix2))
    error('There is at least one non-zero value on one of the diagonals!');
end

% user message if verbose
if verbose
    disp([char(10), 'Called deltaCon function with input args: ',...
        char(10), 'Adjacency (connectivity) matrices of size ', num2str(size(adjMatrix1)),...
        char(10), 'Epsilon (small constant for neighbour weighting): ', num2str(epsilon),...
        char(10), 'Verbosity: ', num2str(verbose), char(10)]);
end


%% Calculate DeltaCon (DeltaCon0 in the papers)

% number of nodes
nodeNo = size(adjMatrix1, 1);

% identity
I = eye(nodeNo);

% node degrees
D1 = diag(sum(adjMatrix1, 1));
D2 = diag(sum(adjMatrix2, 1));

% node affinity matrices
S1 = inv((I+epsilon^2*D1-epsilon*adjMatrix1));
S2 = inv((I+epsilon^2*D2-epsilon*adjMatrix2));

% DeltaCon distance (root Euclidean a.k.a. Matusita distance)
deltaConDist = (sum(sum((S1.^0.5-S2.^0.5).^2, 1), 2))^0.5;

% similarity is bounded to [0 1]
deltaConSim = 1/(1+deltaConDist);

% user message if verbose
if verbose
    disp([char(10), 'DeltaCon distance: ', num2str(deltaConDist),...
        char(10), 'DeltaCon similarity: ', num2str(deltaConSim), char(10)]);
end


return
















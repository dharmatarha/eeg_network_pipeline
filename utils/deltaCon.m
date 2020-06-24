function [deltaConDist, S1, S2] = deltaCon(adjMatrix1, adjMatrix2, varargin)
%% Calculate DeltaCon network similarity measure 
%
% USAGE: [deltaConDist, S1, S2] = deltaCon(adjMatrix1, adjMatrix2, epsilon = 0.01, verbose = true)
% 
% The function calculates the DeltaCon network similarity measure for the
% two networks defined by the adjancency matrices "adjMatrix1" and
% "adjMatrix2". See the details in:
%   Koutra et al. (2013. Deltacon: A principled massive-graph 
%   similarity function. 
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
%               network.
% adjMatrix2    - Numeric square matrix, adjacency matrix for the second 
%               network, same size as "adjMatrix1"
%
% Optional input:
% epsilon       - Numeric value, small constant. Controls the relative 
%               weight of node ranks and the adjacency matrix. See the 
%               paper and the equations for details. Defaults to 0.01.
% verbose       - Logical value (true or false). Verbosity, "false" 
%               meaning no user messages, "true" meaning user messages.
%
% Outputs:
% deltaConDist  - Numeric value, DeltaCon distance. 
% S1            - Numeric matrix, pairwise node affinities for
%               "adjMatrix1". Same size as "adjMatrix1" (nodes X nodes).
% S2            - Numeric matrix, pairwise node affinities for
%               "adjMatrix2". Same size as "adjMatrix2" (nodes X nodes).
%


%% Input checks

% check number of arguments
if ~ismember(nargin, 2:4)
    error(['Function deltaCon requires input args "adMatrix1" and "adMatrix2", ',...
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
        elseif islogical(varargin{v}) && ~exist('verbose', 'var')
            verbose = varargin{v};
        else
            error('An optional input arg does not match nicely to "epsilon" or "verbose"!');
        end
    end
else
    epsilon = 0.01;
    verbose = true;
end









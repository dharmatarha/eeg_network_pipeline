function [withinConn, betweenConn, withinEdgeIndices, betweenEdgeIndices] = modConnections(modules, connMatrix)
%% Calculate within- and between-module connectivities
%
% USAGE: [withinConn, 
%         betweenConn, 
%         withinEdgeIndices, 
%         betweenEdgeIndices] = modConnections(modules, connMatrix)
%
% Calculates within- and between-module connectivity values for a given 
% network (defined by the connectivity / adjacency matrix) and partition
% (as defined by the module membership vector). Between-module connectivity
% is not divided into specific module-pairings, only an aggregate measure
% is returned. 
% The function also returns edge indices (linear indices for the supplied
% connectivity / adjacency matrix) for within- and between module edges. 
% Indices for within-module edges are returned separately for each module.
%
% For UNDIRECTED NETWORKS!
% We only work with upper-triangle connectivity data (above the main 
% diagonal), lower triangle values are ignored.
%
% Inputs:
% modules       - Numeric vector, containing the partition (module
%               membership per node). Might only contain positive integers
%               or zeros. Its length must correspond to the size of
%               "connMatrix" ( length(modules)==size(connmatrix, 1) ).
% connMatrix    - Numeric square matrix, connectivity / adjacency matrix of
%               the network. Only the upper triangle above the main
%               diagonal is taken into account for calculations.
%
% Outputs:
% withinConn         - Numeric vector sized (no. of modules X 1). Its 
%               values are the sums of within-module connectivity values, 
%               one for each module, in ascending order of module indices. 
% betweenConn        - Numeric value, sum of between-module connectivity
%               values.
% withinEdgeIndices  - Cell array sized (no. of modules X 1). Each cell
%               holds a column vector with the linear indices of all 
%               within-module edges. In ascending order of module indices.
% betweenEdgeIndices - Numeric column vector holding the linear indices of
%               all between-module edges.
%


%% Input checks

% no. of args
if nargin ~= 2
    error('Function "modConnections" requires input args "modules" and "connMatrix"!');
end
% arg values
if ~isnumeric(modules) || ~isvector(modules) || any(mod(modules, 1)~=0)
    error('Input arg "modules" should be a numeric vector containing only positive integer and zero values!');
end
if ~isnumeric(connMatrix) || ~ismatrix(connMatrix) || ~isequal(size(connMatrix, 1), size(connMatrix, 2))
    error('Input arg "edges" should be a numeric square matrix!');
end
if length(modules) ~= size(connMatrix, 1)
    error('Input arg "modules" should have equal length to the size of the dimensions of arg "connMatrix"!');
end
% might not even need this...
if ~iscolumn(modules)
    modules = modules';
end


%% Basics

% number of nodes
nodeNo = length(modules);
% get unique module values
modIndices = unique(modules);
% number of modules
moduleNo = length(modIndices);

% preallocate cell array holding edge indices for each module
withinEdgeIndices = cell(moduleNo, 1);
% preallocate vector of within-module connectivity values
withinConn = nan(moduleNo, 1);
% define var holding indices for all within-module edges, across all
% modules
allWithinEdges = [];


%% Looop through modules, extract within-module connections

for m = 1:moduleNo
    
    % get indice for current module
    currentModule = modIndices(m);
    % find corresponding nodes
    nodes = find(modules==currentModule);
    
    % if it even makes sense to talk about within-module connections
    if length(nodes) > 1
        
        % preallocate for var holding end nodes for each within-module edge
        nodePairs = nan(length(nodes)^2, 2);
        
        % loop through all nodes in module
        for n = 1:length(nodes)
            % first we collect all end nodes for all possible edges
            nodePairs((n-1)*length(nodes)+1:n*length(nodes), :) = [repmat(nodes(n), [length(nodes), 1]), nodes];
        end
        
        % then we trim edges not in the upper triangle above the main
        % diagonal
        nodePairs(nodePairs(:, 1) >= nodePairs(:, 2), :) = [];
        % get linear indices for the edges defined by nodePairs
        withinEdgeIndices{m} = sub2ind([nodeNo nodeNo], nodePairs(:, 1), nodePairs(:,2));
        
        % sum of connMatrix values for selected indices = within-module
        % connectivity
        withinConn(m) = sum(connMatrix(withinEdgeIndices{m}));
        
        % aggregate all within-module edge indices into a var
        allWithinEdges = [allWithinEdges; withinEdgeIndices{m}];
        
    else  % if there is only one node in module
        
        withinEdgeIndices{m} = [];  % empty array for within-module edges  
        
    end
    
end
            
            
%% Define across-module connections
% by definition, all edges that are not within-module edges, are 
% between-module ones            

% linear indices for all edges in upper triangle
allEdges = find(triu(true(nodeNo), 1));
% get indices of all between-module edge
betweenEdgeIndices = setdiff(allEdges, allWithinEdges);

% between-module connectivity
betweenConn = sum(connMatrix(betweenEdgeIndices));
            
            


return
            
            
            
            
            
            
            
            
            
            
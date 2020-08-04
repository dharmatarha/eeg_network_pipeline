function nodalRandConn = nodalRandomization(multiLayerConn, omega, nodeNo) 
%% Nodal randomization of mult-layer connectivity matrix
%
% USAGE: nodalRandConn = nodalRandomization(multiLayerConn, omega)
%
% For UNDIRECTED weighted or binary connectivity matrices where the only
% inter-layer connections are between the same nodes in subsequent layers
% (that is, between node "i" in layer "l" and node "i" in layer "l+1").
% Also, we assume that all inter-layer connections are uniform with value
% "omega". This is the usual case for temporally ordered connectivity
% matrices from EEG/MEG/fMRI data.
% 
% The function randomizes inter-layer connections so that an edge from node
% "i" in layer "l" will connect to a random node "j~=i" in layer "l+1".
% This method is in line with the nodel randomization method described in
% Bassett et al., 2013: Robust detection of dynamic comnmunity structure in
% networks
%
% Mandatory inputs:
% multiLayerConn    - Numeric matrix, might be sparse. Sized (node no. X
%                     layer no., node no. X layer no.). Connectivity matrix 
%                     including inter-layer connections, for example as in
%                     the output of getMultiLayerConnMatrix. 
% omega             - Numeric value. Connectivity value of all inter-layer
%                     edges. Could be deduced from input matrix, but is
%                     faster / simpler as an input.
% nodeNo            - Numeric value. Number of nodes in one layer of the 
%                     multi-layer connectivity matrix. As with "omega", 
%                     could be deduced from input matrix, but is faster 
%                     / simpler as an input.
%
% Outputs:
% nodalRandConn     - Numeric matrix, sparse if "multiLayerConn" is sparse.
%                     Sized (node no. X layer no., node no. X layer no.). 
%                     Connectivity matrix with randomized inter-layer 
%                     connections.
%


%% Input checks

if nargin ~= 3
    error('Function nodalRandomization requires input args "multiLayerConn", "omega", "nodeNo"!');
end
if ~isnumeric(multiLayerConn) || ~ismatrix(multiLayerConn) || size(multiLayerConn, 1) ~= size(multiLayerConn, 2)
    error('Input arg "multiLayerConn" should be a square numeric matrix!');
end
if ~isnumeric(omega) || numel(omega)~=1 
    error('Input arg "omega" should be a numeric value!');
end
if ~isnumeric(nodeNo) || numel(nodeNo)~=1 || mod(size(multiLayerConn, 1), nodeNo)~=0
    error('Input arg "nodeNo" should be a numeric value and a divisor of size(multiLayer, 1)!');
end

% get number of layers
layerNo = size(multiLayerConn, 1)/nodeNo;

% loop through layers
for l=1:layerNo-1
    
    % get indices for existing inter-layer connections between layer "l" and "l+1" 
    idxSource = [(l-1)*nodeNo+1:l*nodeNo]';
    idxTarget = idxSource+nodeNo;
    
    % get random target connections
    idxTargetPerm = randi([idxTarget(1), idxTarget(nodeNo)], [nodeNo, 1]);
    % permute until there is none on diagonal
    if any(idxTarget == idxTargetPerm)
        tmpIdx = find(idxTarget==idxTargetPerm);
        for z = tmpIdx'
            tmpTarget = idxTarget;
            tmpTarget(tmpTarget==idxSource(z)+nodeNo) = [];
            idxTargetPerm(z) = tmpTarget(randperm(nodeNo-1, 1));
        end
    end
    
    % create sub-matrix for inter-layer connections
    subM = zeros(nodeNo);
    % add omega values for random connections
    linIdx = sub2ind([nodeNo, nodeNo], idxSource-(l-1)*nodeNo, idxTargetPerm-l*nodeNo);
    subM(linIdx) = omega;
   
    % adjust the sparse connection matrix with the new sub-matrices (subM
    % and its inverse)
    multiLayerConn(idxSource(1):idxSource(nodeNo), idxTarget(1):idxTarget(nodeNo)) = subM;
    multiLayerConn(idxTarget(1):idxTarget(nodeNo), idxSource(1):idxSource(nodeNo)) = subM';
    
end  % for l

% set output var
nodalRandConn = multiLayerConn;


return
    
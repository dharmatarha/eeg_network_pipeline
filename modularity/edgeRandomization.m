function edgeRandConn = edgeRandomization(connectivityMatrix, existingEdges)
%% Performing edge randomization for a matrix
%
% USAGE: edgeRandConn = edgeRandomization(connectivityMatrix, existingEdges=true)
%
% Calculates the randomized version of a given upper triangular 
% connectivity matrix (excluding the main diagonal).
% Calculation is carried out by randomly reassigning of edge weights.
% 
% Mandatory inputs:
% connectivityMatrix      - Numeric matrix. Connectivity matrix for  
 %                          which the randomization must be carried out.  
%                           May contain NaN entries.
%
% Optional inputs:
% existingEdges           - Logical value. Flag for permuting only among
%                           existing edges (that is, 0 or NaN values in  
%                           upper triangle are kept at the same location 
%                           for randomized matrix). Defaults to true.
%
% Output:
% edgeRandConn            - Numeric matrix. Randomized connectivity matrix 
%                           having the same dimensions as the original 
%                           connectivity matrix (arg "connectivityMatrix").
% 


%% Input checks

% Check no. of args
if ~ismembertol(nargin, 1:2)
    error(['Function edgeRandomization requires input arg ',...
        '"connectivityMatrix" while input arg "existingEdges" is optional!']);
end
% check mandatory args
if ~ismatrix(connectivityMatrix)  || size(connectivityMatrix, 1) ~= size(connectivityMatrix, 2)
    error('Input must be a square matrix!');
end
% check optional args 
if nargin == 1
    existingEdges = true;
else
    if ~islogical(existingEdges) || ~isequal(numel(existingEdges), 1)
        error('Optional arg "existingEdges" should be a logical value!');
    end
end


%% Edge randomization

% get channel no.
[numberOfChannels, ~] = size(connectivityMatrix);

% some value that should never occur as a connectivity value in our use
% cases
blackSwan = -999;
% complain if - for whatever unimaginable reason - any matrix element is exactly our black swan 
if any(any(connectivityMatrix==blackSwan))
    error('The input matrix contains a value we thought we would never encounter and use internally in the script...');
end

% Set lower triangular + main diagonal of the connectivity matrix to some
% impossible value
connectivityMatrix(tril(true(numberOfChannels))) = blackSwan;  

% if we only permute among existing edges, get a mask for nonzero & non-NaN
% value in upper triangle
if existingEdges
    mask = triu(true(numberOfChannels), 1);
    mask(isnan(connectivityMatrix)) = false;
    mask(connectivityMatrix == 0) = false;
end

% extract connectivity values from upper triangle into a vector
connValuesLin = connectivityMatrix(:);  % linearize values of the matrix
connValuesLin(connValuesLin==blackSwan) = [];  % delete blackSwan values - supposed to be only the lower triangle + diagonal
% delete also zero & NaN values if permutation is only among existing edges
if existingEdges
    connValuesLin(isnan(connValuesLin)) = [];
    connValuesLin(connValuesLin == 0) = [];
end

% randomize the extracted connectivity values
randConnValues = connValuesLin(randperm(length(connValuesLin)));

% fill an upper triangle of a matrix (numberOfChannels,
% numberOfChannels) with the randomized values
edgeRandConn = nan(numberOfChannels);
% if permutation was only among existing edges, assignment is via the mask
if existingEdges
    edgeRandConn(mask) = randConnValues;
    % switch upper triangle NaNs to zeros if there were any in original
    % connectivity matrix
    edgeRandConn(connectivityMatrix==0) = 0;
else
    edgeRandConn(triu(true(numberOfChannels), 1)) = randConnValues;
end



return


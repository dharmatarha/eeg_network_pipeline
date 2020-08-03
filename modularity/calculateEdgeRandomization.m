function [randomizedConnectivityMatrix] = calculateEdgeRandomization(connectivityMatrix)
%% Performing edge randomization for a matrix
%
% USAGE: randomizedConnectivityMatrix = calculateEdgeRandomization(connectivityMatrix)
%
% Calculates the randomized version of a given upper triangular connectivity matrix (excluding the main diagonal).
% Calculation is carried out by random reassigning of edge weights.
% 
% Mandatory inputs:
% connectivityMatrix            - Connectivity matrix for which the randomization must be carried out.
%                                 Must be a matrix. May contain NaN entries.
%
% Optional inputs:
% None.
%
% Output:
% randomizedConnectivityMatrix - Randomized connectivity matrix having the same dimensions as the original connectivity matrix.
% 


%% Input checks

% Check for mandatory arguments
if nargin ~= 1
    error('Measured connectivity matrix is required!');
end
if ~ismatrix(connectivityMatrix)
    error('Input must be a matrix!');
end
if size(connectivityMatrix, 1) ~= size(connectivityMatrix, 2)
    error('Input must be a square matrix!');
end


%% Vectorized version:

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

% extract connectivity values from upper triangle into a vector
connValuesLin = connectivityMatrix(:);  % linearize values of the matrix
connValuesLin(connValuesLin==blackSwan) = [];  % delete blackSwan values - supposed to be only the lower triangle + diagonal

% randomize the extracted connectivity values
randConnValues = connValuesLin(randperm(length(connValuesLin)));

% fill an upper triangle of a matrix (numberOfChannels,
% numberOfChannels) with the randomized values
randomizedConnectivityMatrix = nan(numberOfChannels);
randomizedConnectivityMatrix(triu(true(numberOfChannels), 1)) = randConnValues;

return

% %% Determine data dimensions  
% 
% [numberOfChannels, ~] = size(connectivityMatrix);
% % preallocate the randomized matrix
% randomizedConnectivityMatrix = zeros(numberOfChannels, numberOfChannels);
% % Set to lower triangular part and to main diagonal of the connectivity matrix to NaN
% for nanColumnIndex = 1 : numberOfChannels
%     for nanRowIdex = nanColumnIndex : numberOfChannels
%         connectivityMatrix(nanRowIdex, nanColumnIndex) = NaN;
%     end
% end
% numberOfNaNElements = sum(sum(isnan(connectivityMatrix)));
% numberOfEdgeWeightsToReassign = numel(connectivityMatrix) - numberOfNaNElements;
% 
% %% Calculate the randomized matrix
% 
% numberOfEdgeWeightsAlreadyReassigned = 0;
% indicesOfEdgesAlreadyReassigned = zeros(1, numberOfEdgeWeightsToReassign);
% if mod(numberOfEdgeWeightsToReassign, 2) == 1
%     numberOfEdgeWeightsToReassign = numberOfEdgeWeightsToReassign - 1;
% end
% while numberOfEdgeWeightsAlreadyReassigned < numberOfEdgeWeightsToReassign
%     firstEdgeFirstNodeIndex = randi([1 numberOfChannels-1]);
%     firstEdgeSecondNodeIndex = randi([firstEdgeFirstNodeIndex+1 numberOfChannels]); % we consider only upper triangular part of the matrix
%     secondEdgeFirstNodeIndex = randi([1 numberOfChannels-1]);
%     secondEdgeSecondNodeIndex = randi([secondEdgeFirstNodeIndex+1 numberOfChannels]); % we consider only upper triangular part of the matrix
%     
%     firstEdgeIndex = (firstEdgeFirstNodeIndex-1) * numberOfChannels + firstEdgeSecondNodeIndex;
%     secondEdgeIndex = (secondEdgeFirstNodeIndex-1) * numberOfChannels + secondEdgeSecondNodeIndex;
%     if firstEdgeIndex ~= secondEdgeIndex
%         if (~isnan(connectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex))) && (~isnan(connectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex)))
%             if (~ismember(firstEdgeIndex, indicesOfEdgesAlreadyReassigned)) && (~ismember(secondEdgeIndex, indicesOfEdgesAlreadyReassigned))
%                 randomizedConnectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex) = connectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex);
%                 randomizedConnectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex) = connectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex);
%                 indicesOfEdgesAlreadyReassigned(numberOfEdgeWeightsAlreadyReassigned+1) = firstEdgeIndex;
%                 indicesOfEdgesAlreadyReassigned(numberOfEdgeWeightsAlreadyReassigned+2) = secondEdgeIndex;
%                 numberOfEdgeWeightsAlreadyReassigned = numberOfEdgeWeightsAlreadyReassigned + 2;
%             end
%         end
%     end
% end
% 
% NaNElements = isnan(connectivityMatrix);
% randomizedConnectivityMatrix(NaNElements) = NaN;
% edgeWeightsLeftOut = find(randomizedConnectivityMatrix == 0);
% randomizedConnectivityMatrix(edgeWeightsLeftOut) = connectivityMatrix(edgeWeightsLeftOut);
% 
% end
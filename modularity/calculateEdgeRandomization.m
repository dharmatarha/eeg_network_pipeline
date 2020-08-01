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
%%

%% Input checks

% Check for mandatory arguments
if nargin < 1
    error('Measured connectivity matrix is required!');
end
if (numel(size(connectivityMatrix)) ~= 2)
    error('Input must be a matrix!');
end
if size(connectivityMatrix, 1) ~= size(connectivityMatrix, 2)
    error('Input must be a matrix!');
end


%% Determine data dimensions  

[numberOfChannels, ~] = size(connectivityMatrix);
% preallocate the randomized matrix
randomizedConnectivityMatrix = zeros(numberOfChannels, numberOfChannels);
% Set to lower triangular part and to main diagonal of the connectivity matrix to NaN
for nanColumnIndex = 1 : numberOfChannels
    for nanRowIdex = nanColumnIndex : numberOfChannels
        connectivityMatrix(nanRowIdex, nanColumnIndex) = NaN;
    end
end
numberOfNaNElements = sum(sum(isnan(connectivityMatrix)));
numberOfEdgeWeightsToReassign = numel(connectivityMatrix) - numberOfNaNElements;

%% Calculate the randomized matrix

numberOfEdgeWeightsAlreadyReassigned = 0;
indicesOfEdgesAlreadyReassigned = zeros(1, numberOfEdgeWeightsToReassign);
if mod(numberOfEdgeWeightsToReassign, 2) == 1
    numberOfEdgeWeightsToReassign = numberOfEdgeWeightsToReassign - 1;
end
while numberOfEdgeWeightsAlreadyReassigned < numberOfEdgeWeightsToReassign
    firstEdgeFirstNodeIndex = randi([1 numberOfChannels-1]);
    firstEdgeSecondNodeIndex = randi([firstEdgeFirstNodeIndex+1 numberOfChannels]); % we consider only upper triangular part of the matrix
    secondEdgeFirstNodeIndex = randi([1 numberOfChannels-1]);
    secondEdgeSecondNodeIndex = randi([secondEdgeFirstNodeIndex+1 numberOfChannels]); % we consider only upper triangular part of the matrix
    
    firstEdgeIndex = (firstEdgeFirstNodeIndex-1) * numberOfChannels + firstEdgeSecondNodeIndex;
    secondEdgeIndex = (secondEdgeFirstNodeIndex-1) * numberOfChannels + secondEdgeSecondNodeIndex;
    if firstEdgeIndex ~= secondEdgeIndex
        if (~isnan(connectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex))) && (~isnan(connectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex)))
            if (~ismember(firstEdgeIndex, indicesOfEdgesAlreadyReassigned)) && (~ismember(secondEdgeIndex, indicesOfEdgesAlreadyReassigned))
                randomizedConnectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex) = connectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex);
                randomizedConnectivityMatrix(secondEdgeFirstNodeIndex, secondEdgeSecondNodeIndex) = connectivityMatrix(firstEdgeFirstNodeIndex, firstEdgeSecondNodeIndex);
                indicesOfEdgesAlreadyReassigned(numberOfEdgeWeightsAlreadyReassigned+1) = firstEdgeIndex;
                indicesOfEdgesAlreadyReassigned(numberOfEdgeWeightsAlreadyReassigned+2) = secondEdgeIndex;
                numberOfEdgeWeightsAlreadyReassigned = numberOfEdgeWeightsAlreadyReassigned + 2;
            end
        end
    end
end

NaNElements = isnan(connectivityMatrix);
randomizedConnectivityMatrix(NaNElements) = NaN;
edgeWeightsLeftOut = find(randomizedConnectivityMatrix == 0);
randomizedConnectivityMatrix(edgeWeightsLeftOut) = connectivityMatrix(edgeWeightsLeftOut);

end
function [postOptimizationNullNetwork] = calculatePostOptimizationNullNetwork(prunedConnectivity, varargin)
%% Multilayer post-optimization null network calculation
%
% USAGE: postOptimizationNullNetwork = calculatePostOptimizationNullNetwork(prunedConnectivity, intraLayerEdgeWeightRandomization = 1, temporalRandomization = 0)
%
% Calculates the post-optimization null network of connectivity data consisting of epochs and sessions (stories).
% Calculation is carried out using given measured connectivity data by intra-layer edge weight randomization (random reassigning of edge weights) 
% and temporal randomization (random reshuffling of network layers).
% 
% Mandatory inputs:
% prunedConnectivity                - Measured connectivity data containing only significant connections.
%                                     Must have the dimensions (numberOfChannels, numberOfChannels, numberOfEpochs, numberOfStories).
%
% Optional inputs:
% intraLayerEdgeWeightRandomization - Variable indicating if intra-layer edge weight randomization (random reassigning of edge weights) should be performed.
%                                     If empty, default value = 1 is used.
% temporalRandomization             - Variable indicating if temporal randomization (random reshuffling of network layers) should be performed.
%                                     If empty, default value = 0 is used.
%
% Output:
% postOptimizationNullNetwork       - Post-optimization null network having the same dimensions as the original network (numberOfChannels, numberOfChannels, numberOfEpochs, numberOfStories).
% 
%%

%% Input checks

% Check for mandatory arguments
if nargin < 1
    error('Measured connectivity data is required!');
end
if (numel(size(prunedConnectivity)) ~= 4)
    error('Measured connectivity data dimensions are not as required!');
end

% Check optional arguments
if ~isempty(varargin)
    % If too many, raise error
    if length(varargin) > 2
        error('Too many variable inputs. Only "intraLayerEdgeWeightRandomization" and "temporalRandomization" are allowed!');
    % If only one, it is intraLayerEdgeWeightRandomization
    elseif length(varargin) == 1
            intraLayerEdgeWeightRandomization = varargin{1};
            temporalRandomization = 0;
    % If two, the first one is intraLayerEdgeWeightRandomization, the second one is temporalRandomization
    elseif length(varargin) == 2
        intraLayerEdgeWeightRandomization = varargin{1};
        temporalRandomization = varargin{2};
    end
else
    intraLayerEdgeWeightRandomization = 1;
    temporalRandomization = 0;
end


%% Determine data dimensions  

[numberOfChannels, ~, numberOfEpochs, numberOfStories] = size(prunedConnectivity);
% preallocate the tensor for the post-optimization null network
postOptimizationNullNetwork = zeros(numberOfChannels, numberOfChannels, numberOfEpochs, numberOfStories);


%% Calculate one post-optimization null network for each story

for storyIndex = 1 : numberOfStories
    
    if intraLayerEdgeWeightRandomization == 1
        for epochIndex = 1 : numberOfEpochs
            connectivityMatrix = squeeze(prunedConnectivity(:, :, epochIndex, storyIndex));
            randomizedConnectivityMatrix = calculateEdgeRandomization(connectivityMatrix);
            postOptimizationNullNetwork(:, :, epochIndex, storyIndex) = randomizedConnectivityMatrix;       
        end
    else
        postOptimizationNullNetwork(:, :, :, storyIndex) = prunedConnectivity(:, :, :, storyIndex);
    end
    
    if temporalRandomization == 1
        postOptimizationNullNetwork(:, :, :, storyIndex) = calculateTemporalRandomization(squeeze(postOptimizationNullNetwork(:, :, :, storyIndex)));
    end 
    
end

end
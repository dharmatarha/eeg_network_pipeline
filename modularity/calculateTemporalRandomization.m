function [temporallyRandomizedMultiLayerNetwork] = calculateTemporalRandomization(multiLayerNetwork)
%% Performing temporal randomization of a multilayer network
%
% USAGE: temporallyRandomizedMultiLayerNetwork = calculateTemporalRandomization(multiLayerNetwork)
%
% Calculates the temporally randomized version of a given multilayer network.
% Calculation is carried out by random reshuffling of network layers.
% 
% Mandatory inputs:
% multiLayerNetwork                     - Multilayer network for which the randomization must be carried out.
%                                         Must be a 3D tensor. Must have the dimensions (numberOfChannels, numberOfChannels, numberOfLayers).
%
% Optional inputs:
% None.
%
% Output:
% temporallyRandomizedMultiLayerNetwork - Temporally randomized multilayer network having the same dimensions as the original multilayer network.
% 
%

%% Input checks

% Check for mandatory arguments
if nargin < 1
    error('Multilayer network is required!');
end
if (numel(size(multiLayerNetwork)) ~= 3)
    error('Input must be a 3D tensor!');
end

%% Randomize layers

[~, ~, numberOfLayers] = size(multiLayerNetwork);
temporallyRandomizedMultiLayerNetwork = multiLayerNetwork(:, :, randperm(numberOfLayers));

return

% %% Determine data dimensions  
% 
% [numberOfChannels, ~, numberOfLayers] = size(multiLayerNetwork);
% % preallocate the randomized tensor
% temporallyRandomizedMultiLayerNetwork = zeros(numberOfChannels, numberOfChannels, numberOfLayers);
% 
% 
% %% Calculate the temporally randomized multilayer network
% 
% randomizedLayerIndices = randperm(numberOfLayers);
% 
% for indexOfIndices = 1 : numberOfLayers
%     temporallyRandomizedMultiLayerNetwork(:, :, indexOfIndices) = multiLayerNetwork(:, :, randomizedLayerIndices(indexOfIndices));
% end
% 
% end
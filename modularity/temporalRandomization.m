function tempRandConn = temporalRandomization(multiLayerConn)
%% Performing temporal randomization of a multilayer network
%
% USAGE: tempRandConn = temporalRandomization(multiLayerConn)
%
% Calculates the temporally randomized version of a given multilayer network.
% Calculation is carried out by random reshuffling of network layers.
% 
% Mandatory inputs:
% multiLayerConn     - Numeric tensor, 3D. Multi-layer connectivity
%                      network to be randomized. We assume that
%                      subsequent layers represent ordinal (temporal)
%                      order. Must have dimensions 
%                      (numberOfChannels, numberOfChannels, numberOfLayers).
%
% Optional inputs:
% None.
%
% Output:
% tempRandConn      - Temporally randomized multi-layer connectivity 
%                     network having the same dimensions as the original 
%                     network (arg "multiLayerConn").
% 


%% Input checks

% Check for mandatory arguments
if nargin ~= 1
    error('Input arg "multiLayerConn" is required!');
end
if (numel(size(multiLayerConn)) ~= 3)
    error('Input must be a 3D tensor!');
end


%% Randomize layers

[~, ~, numberOfLayers] = size(multiLayerConn);
tempRandConn = multiLayerConn(:, :, randperm(numberOfLayers));


return


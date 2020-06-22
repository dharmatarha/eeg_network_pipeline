function [B] = calculateMultiLayerModularityMatrix(prunedConnectivity, meanSurrogateData, varargin)
%% Multilayer modularity matrix calculation
%
% USAGE: B = calculateMultiLayerModularityMatrix(prunedConnectivity, meanSurrogateData, gamma = 1, omega = 1)
%
% Calculates the multilayer modularity matrices of connectivity data consisting of epochs and sessions (stories).
% Calculation is carried out using given measured connectivity data, surrogate connectivity data, resolution parameter (gamma) and inter-layer edge weight (omega).
% 
% Mandatory inputs:
% prunedConnectivity  - Measured connectivity data containing only significant connections.
%                     Must have the dimensions (numberOfChannels, numberOfChannels, numberOfEpochs, numberOfStories).
% meanSurrogateData   - Mean surrogate connectivity data.
%                     Must have the dimensions (numberOfChannels, numberOfChannels, numberOfEpochs, numberOfStories).
%
% Optional inputs:
% gamma               - Resolution parameter.
%                     If empty, default value = 1 is used.
% omega               - Inter-layer edge weight. Uniform inter-layer connections are assumed (all inter-layer edges have the same weight).
%                     If empty, default value = 1 is used.
%
% Output:
% B                   - Cell array of multilayer modularity matrices.
%                     The number of matrices in the array is equal to the number of sessions (stories).
% 
%%

%% Input checks

% Check for mandatory arguments
if nargin < 2
    error('Measured and surrogate connectivity data is required!');
end
if (numel(size(prunedConnectivity)) ~= 4) || (numel(size(meanSurrogateData)) ~= 4)
    error('Measured or surrogate connectivity data dimensions are not as required!');
end

% Check optional arguments
if ~isempty(varargin)
    % If too many, raise error
    if length(varargin) > 2
        error('Too many variable inputs. Only "gamma" and "omega" are allowed!');
    % If only one, it is gamma
    elseif length(varargin) == 1
            gamma = varargin{1};
            omega = 1;
    % If two, the first one is gamma, the second one is omega
    elseif length(varargin) == 2
        gamma = varargin{1};
        omega = varargin{2};
    end
else
    gamma = 1;
    omega = 1;
end


%% Determine data dimensions  

[numberOfChannels, ~, numberOfEpochs, numberOfStories] = size(prunedConnectivity);
% preallocate the cell array holding one multi-layer matrix per story /
% stimulus
B = cell(numberOfStories, 1);


%% Calculate one multilayer modularity matrix for each story
for storyIndex = 1 : numberOfStories
    % Tensor for intra-layer modularity matrices
    intraLayerModularityTensor = cell(numberOfEpochs, numberOfChannels, numberOfChannels);
    
    for epochIndex = 1 : numberOfEpochs
        % Measured and surrogate connectivity matrices for a given epoch
        connectivityMatrix = squeeze(prunedConnectivity(:, :, epochIndex, storyIndex));
        surrogateMatrix = squeeze(meanSurrogateData(:, :, epochIndex, storyIndex));
        
        % Keep only connections in both measured and surrogate matrices, where measured connectivity exists
        surrogateMatrix(isnan(connectivityMatrix)) = 0;
        connectivityMatrix(isnan(connectivityMatrix)) = 0;
        
        % Normalize measured and surrogate connectivity matrices
        normalizedConnectivityMatrix = normalizeMatrix(connectivityMatrix);
        normalizedSurrogateMatrix = normalizeMatrix(surrogateMatrix);
        % Calculate intra-layer modularity matrix
        intraLayerModularityMatrix = normalizedConnectivityMatrix - gamma * normalizedSurrogateMatrix;
        
        % If the intra-layer modularity matrix is not symmetric, symmetrization is forced
        if nnz(intraLayerModularityMatrix-intraLayerModularityMatrix')
            intraLayerModularityMatrix = (intraLayerModularityMatrix + intraLayerModularityMatrix')/2; %disp('WARNING: Forced symmetric intra-layer modularity matrix ')
        end
        
        intraLayerModularityTensor{epochIndex} = intraLayerModularityMatrix;
    end
    
    N = numberOfChannels;
    T = numberOfEpochs;
    % Allocate sparse matrix for the multilayer modularity matrix of the given story
    B_givenStory = spalloc(N*T, N*T, N*N*T+2*N*T);
    % Iterate over all layers
    for s = 1:T
        indx = [1:N] + (s-1)*N;
        % Diagonal blocks are the intra-layer modularity matrices
        B_givenStory(indx, indx) = intraLayerModularityTensor{s};
    end
    % Off-diagonal blocks contain inter-layer connections
    B_givenStory = B_givenStory + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
    
    B{storyIndex} = B_givenStory;
end

end
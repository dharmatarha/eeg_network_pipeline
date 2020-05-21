function [B] = calculateMultiLayerModularityMatrix(prunedConnectivity, meanSurrogateData, storyIndex, gamma, omega)

% numberOfStories = size(prunedConnectivity, 4);
numberOfEpochs = size(prunedConnectivity, 3);
numberOfChannels = size(prunedConnectivity, 1);

multiLayerModularityTensor = cell(numberOfEpochs, numberOfChannels, numberOfChannels);

for epochIndex = 1 : numberOfEpochs
    connectivityMatrix = squeeze(prunedConnectivity(:, :, epochIndex, storyIndex));
    surrogateMatrix = squeeze(meanSurrogateData(:, :, epochIndex, storyIndex));
    
    for elementIndex = 1 : numel(connectivityMatrix)
        if isnan(connectivityMatrix(elementIndex))
            connectivityMatrix(elementIndex) = 0;
            surrogateMatrix(elementIndex) = 0;
        end
    end
    
    normalizedConnectivityMatrix = normalizeMatrix(connectivityMatrix);
    normalizedSurrogateMatrix = normalizeMatrix(surrogateMatrix);
    intraLayerModularityMatrix = normalizedConnectivityMatrix - gamma * normalizedSurrogateMatrix;
    
    if nnz(intraLayerModularityMatrix-intraLayerModularityMatrix')
        intraLayerModularityMatrix = (intraLayerModularityMatrix + intraLayerModularityMatrix')/2; disp('WARNING: Forced symmetric intra-layer modularity matrix ')
    end
    
    multiLayerModularityTensor{epochIndex} = intraLayerModularityMatrix;
end

N = numberOfChannels;
T = numberOfEpochs;
B = spalloc(N*T, N*T, N*N*T+2*N*T);
for s = 1:T
    indx = [1:N] + (s-1)*N;
    B(indx, indx) = multiLayerModularityTensor{s};
end
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

end
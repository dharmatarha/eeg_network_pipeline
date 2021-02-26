clear all;
close all;

pathName = ('D:\mult_dim_scale\');

for index = 1 : 4
    switch index
        case 1
            fileName = ('group_alpha_plv.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
            connectivityMetricString = 'PLV ';
        case 2
            fileName = ('group_alpha_iplv.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
            connectivityMetricString = 'iPLV ';
        case 3
            fileName = ('group_alpha_ampCorr.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
            connectivityMetricString = 'ampCorr ';
        case 4
            fileName = ('group_alpha_orthAmpCorr.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
            connectivityMetricString = 'orthAmpCorr ';
        otherwise
    end
    
    [distRes_corr] = connDistanceTest_betweenSubject(connectivityTensor, 'corr');
    [distRes_eucl] = connDistanceTest_betweenSubject(connectivityTensor, 'eucl');
    [distRes_adjacencySpectral] = connDistanceTest_betweenSubject(connectivityTensor, 'adjacencySpectral');
    [distRes_LaplacianSpectral] = connDistanceTest_betweenSubject(connectivityTensor, 'LaplacianSpectral');
    [distRes_deltaCon] = connDistanceTest_betweenSubject(connectivityTensor, 'deltaCon');
    
    distRes_corr(isnan(distRes_corr)) = 0;
    distRes_corr = distRes_corr + distRes_corr';
    distRes_eucl(isnan(distRes_eucl)) = 0;
    distRes_eucl = distRes_eucl + distRes_eucl';
    distRes_adjacencySpectral(isnan(distRes_adjacencySpectral)) = 0;
    distRes_adjacencySpectral = distRes_adjacencySpectral + distRes_adjacencySpectral';
    distRes_LaplacianSpectral(isnan(distRes_LaplacianSpectral)) = 0;
    distRes_LaplacianSpectral = distRes_LaplacianSpectral + distRes_LaplacianSpectral';
    distRes_deltaCon(isnan(distRes_deltaCon)) = 0;
    distRes_deltaCon = distRes_deltaCon + distRes_deltaCon';
    
    [Y_corr, eigvals_corr] = cmdscale(distRes_corr);
    [Y_eucl, eigvals_eucl] = cmdscale(distRes_eucl);
    [Y_adjacencySpectral, eigvals_adjacencySpectral] = cmdscale(distRes_adjacencySpectral);
    [Y_LaplacianSpectral, eigvals_LaplacianSpectral] = cmdscale(distRes_LaplacianSpectral);
    [Y_deltaCon, eigvals_deltaCon] = cmdscale(distRes_deltaCon);
    
    cumulativeVariance_corr = zeros(1, numel(eigvals_corr));
    cumulativeVariance_eucl = zeros(1, numel(eigvals_eucl));
    cumulativeVariance_adjacencySpectral = zeros(1, numel(eigvals_adjacencySpectral));
    cumulativeVariance_LaplacianSpectral = zeros(1, numel(eigvals_LaplacianSpectral));
    cumulativeVariance_deltaCon = zeros(1, numel(eigvals_deltaCon));
    for eigenvalueIndex = 1 : numel(eigvals_corr)
        cumulativeVariance_corr(eigenvalueIndex) = sum(eigvals_corr(1:eigenvalueIndex)) / sum(eigvals_corr);
        cumulativeVariance_eucl(eigenvalueIndex) = sum(eigvals_eucl(1:eigenvalueIndex)) / sum(eigvals_eucl);
        cumulativeVariance_adjacencySpectral(eigenvalueIndex) = sum(eigvals_adjacencySpectral(1:eigenvalueIndex)) / sum(eigvals_adjacencySpectral);
        cumulativeVariance_LaplacianSpectral(eigenvalueIndex) = sum(eigvals_LaplacianSpectral(1:eigenvalueIndex)) / sum(eigvals_LaplacianSpectral);
        cumulativeVariance_deltaCon(eigenvalueIndex) = sum(eigvals_deltaCon(1:eigenvalueIndex)) / sum(eigvals_deltaCon);
    end
    
    figure();
    subplot(3, 2, 1);
    plot(eigvals_corr/sum(eigvals_corr), 'LineWidth', 2);
    title([connectivityMetricString 'Correlation']);
    subplot(3, 2, 2);
    plot(eigvals_eucl/sum(eigvals_eucl), 'LineWidth', 2);
    title([connectivityMetricString 'Eucledian distance']);
    subplot(3, 2, 3);
    plot(eigvals_adjacencySpectral/sum(eigvals_adjacencySpectral), 'LineWidth', 2);
    title([connectivityMetricString 'Adjacency spectral distance']);
    subplot(3, 2, 4);
    plot(eigvals_LaplacianSpectral/sum(eigvals_LaplacianSpectral), 'LineWidth', 2);
    title([connectivityMetricString 'Laplacian spectral distance']);
    subplot(3, 2, 5);
    plot(eigvals_deltaCon/sum(eigvals_deltaCon), 'LineWidth', 2);
    title([connectivityMetricString 'DeltaCon']);
    
    figure();
    subplot(3, 2, 1);
    plot(cumulativeVariance_corr, 'LineWidth', 2);
    ylim([0 1]);
    title([connectivityMetricString 'Correlation']);
    subplot(3, 2, 2);
    plot(cumulativeVariance_eucl, 'LineWidth', 2);
    ylim([0 1]);
    title([connectivityMetricString 'Eucledian distance']);
    subplot(3, 2, 3);
    plot(cumulativeVariance_adjacencySpectral, 'LineWidth', 2);
    ylim([0 1]);
    title([connectivityMetricString 'Adjacency spectral distance']);
    subplot(3, 2, 4);
    plot(cumulativeVariance_LaplacianSpectral, 'LineWidth', 2);
    ylim([0 1]);
    title([connectivityMetricString 'Laplacian spectral distance']);
    subplot(3, 2, 5);
    plot(cumulativeVariance_deltaCon, 'LineWidth', 2);
    ylim([0 1]);
    title([connectivityMetricString 'DeltaCon']);
    
end
    


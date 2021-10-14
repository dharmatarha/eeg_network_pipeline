%% Script for Figure 3A and 3B
%
% The script generates 
% (1) The edge contribution matrix of edges 
% (cohen D values on the connectivity map). 
% Reordered for easy-to-read labeling.
% (2) The edge contribution effect size histogram with the cut-off
% determined by the FDR on p values
%
% 

%% Base params

method = 'ciplv'; 
freq = 'alpha';

simMetric = 'corr';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'v2thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load edge contributions data, determined by "method", "freq" and "thr"

edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_groupv2thrSub_onlyPos.mat']);
data = load(edgeContrF);

% sanity check - were the edge contributions calculated with the requested
% similarity metric?
if ~strcmp(simMetric, data.simMetric)
    error(['Bad similarity metric! Edge contributions file specified simMetric: ', data.simMetric]);
end


%% Get cohen d values and p values in connectivity map format

ps = data.permRes.pEstAll; ps(ps==0) = 1 / data.permNo;
ds = data.permRes.cohendAll;

% in connectiviy map formats, by populating the upper triangles
pMap = nan(62); pMap(triu(true(62), 1)) = ps;
dMap = nan(62); dMap(triu(true(62), 1)) = ds;


%% Get cohen d (effect size) map as heatmap - fig 3A


% Load ROI label names
labels = load('/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;  % short names


% Get reordered heatmap
[connMatrixRo, labelsRo, h, hs] = heatmapReordered(dMap, labels);

% Save out Fig 3A
saveas(gcf, fullfile(baseDir, ['edgeEffectSize_heatmap_', thr, '_', freq, '_', method, '.svg']));
saveas(gcf, fullfile(baseDir, ['edgeEffectSize_heatmap_', thr, '_', freq, '_', method, '.png']));
close(gcf);


%% Get histogram of effect sizes with statistical cut-off (FDR)

% fdr
q = 0.01;
fdrMethod = 'bh';
[h, crit] = fdr(ps, q, fdrMethod);

% user message
disp([char(10), char(10),...
    'Number of significant edges: ', num2str(sum(h)), char(10), ...
    'Critical p value: ', num2str(crit), char(10)]);
    
% test if there is a clear cutoff for effect sizes as well
if any(ds(~h) > min(ds(h)))
    histFlag = false;
    tmpNo = sum(ds(~h) > min(ds(h)));
    disp([char(10), 'There is no clean cutoff in terms of effect size.', ...
        char(10), 'There are ', num2str(tmpNo), ' edges with relatively ',...
        'large effect sizes that do not survive the FDR.', char(10)]);
else
    histFlag = true;
    dCutoff = min(ds(h));
    disp([char(10), 'There is a clear cutoff in terms of effect size.',...
        char(10), 'Cutoff (cohen D): ', num2str(dCutoff), char(10)]);
end


% histogram
if histFlag
    
    % plotting params
    gcfMainPos = [0, 0, 0.40, 0.60];
    gcfBackground = [1 1 1];
    fontSize = 14;
    fontWeight = 'Bold';    
    facecolor = [0.7290, 0.5440, 0.0950];
    linewidth = 1.5;
    
    hFig = histogram(ds, 40);
    xlabel('Effect size (Cohen''s d)');
    ylabel('Count');
    hFig.FaceColor = facecolor;
    hFig.LineWidth = linewidth;
    
    % background to white
    set(gcf,'Color', gcfBackground);
    % set size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
    % Font size, weight
    set(gca, 'FontSize', fontSize);
    set(gca, 'FontWeight', fontWeight);
    
    % cutoff line
    l = line([dCutoff, dCutoff], [0 250], 'LineWidth', 2);
    
    saveas(gcf, fullfile(baseDir, ['edgeEffectSize_histogram_', thr, '_', freq, '_', method, '.svg']));
    saveas(gcf, fullfile(baseDir, ['edgeEffectSize_histogram_', thr, '_', freq, '_', method, '.png']));
    close(gcf);

end






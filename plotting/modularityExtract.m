% load connectivity results
resFile = '/media/adamb/bonczData/hyperscan/edgePruningResults_iplv/alpha/avg_alpha_edgePruningInfo.mat';
r = load(resFile);
connData = r.meanConn.prunedConnAvg(:,:,:,1);

% load modularity results
modularityResFile = '/media/adamb/bonczData/hyperscan/modularityResults/tensorNormRes/alpha_stim1_iplv_multiCommWrapperRes_real_iterated_postpr.mat';
m = load(modularityResFile);

% get modularity around target gamma, omega pairing
targets = [1, 1];
[~, gIdx] = min(abs(m.gammaValues-targets(1)));
[~, oIdx] = min(abs(m.omegaValues-targets(2)));

% get modularity
mod = m.res.consSim(gIdx, oIdx, :);
mod = double(squeeze(mod));

% reshape
modAll = reshape(mod, [62, 40]);

% labels
roiFile = '/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat';
roi = load(roiFile);
labels = roi.rois;

% colors
colorFile = '/home/adamb/eeg_network_pipeline/utils/colorTriplets.mat';
c = load(colorFile);
colorTriplets = c.colorTriplets;

% thr
trimmingThr = [0.02, 0];

% rearrange label names and modularity indices
% expected set of labels
newLabels = {'lateralorbitofrontal L', 'medialorbitofrontal L', 'parsorbitalis L', 'parstriangularis L', 'parsopercularis L', 'rostralmiddlefrontal L', 'caudalmiddlefrontal L', 'superiorfrontal L', 'precentral L',...
    'rostralanteriorcingulate L', 'caudalanteriorcingulate L', 'posteriorcingulate L', 'isthmuscingulate L',...
    'transversetemporal L', 'superiortemporal L', 'middletemporal L', 'inferiortemporal L', 'entorhinal L', 'parahippocampal L', 'fusiform L', 'insula L',...
    'supramarginal L', 'inferiorparietal L', 'superiorparietal L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'cuneus L', 'lingual L', 'pericalcarine L', 'lateraloccipital L',...
    'lateraloccipital R', 'pericalcarine R', 'lingual R', 'cuneus R',...
    'precuneus R', 'paracentral R', 'postcentral R', 'superiorparietal R', 'inferiorparietal R', 'supramarginal R',...
    'insula R', 'fusiform R', 'parahippocampal R', 'entorhinal R', 'inferiortemporal R', 'middletemporal R', 'superiortemporal R', 'transversetemporal R',...
    'isthmuscingulate R', 'posteriorcingulate R', 'caudalanteriorcingulate R', 'rostralanteriorcingulate R',...
    'precentral R', 'superiorfrontal R', 'caudalmiddlefrontal R', 'rostralmiddlefrontal R', 'parsopercularis R', 'parstriangularis R', 'parsorbitalis R', 'medialorbitofrontal R', 'lateralorbitofrontal R',...
    };
shiftLabels = 15;
newLabels = [newLabels(end-shiftLabels:end), newLabels(1:end-shiftLabels-1)];

% get intersect of supplied labels and expected ones
c = intersect(labels, newLabels);

% if the two sets match, set data transformation flag
if isequal(c, labels)
    transformFlag = 1;  % transformation flag
else
    error('Unexpected labels');
end


for layerIdx = 1:40
    
    % title
    figTitle = ['Alpha, layer ', num2str(layerIdx), ', stim 1'];
    
    % connectivity
    connMatrix = squeeze(connData(:,:,layerIdx));
    % modularity indices
    modIndicesVector = modAll(:,layerIdx);

    % transform connectivity and module data in case of special label set
    if transformFlag   
    
        % rearrange connectivity matrix based on new ROI/node label
        % order
        [connMatrix, old2new] = matrixReorder(connMatrix, labels, newLabels);
        % apply the same re-ordering to ROI/node module indices
        modIndicesVector = modIndicesVector(old2new);

        % plot
        [mainFig, subFig] = circleGraphPlot(connMatrix, modIndicesVector, colorTriplets, trimmingThr,... 
                                              newLabels, figTitle);
                                          
    end
    
    % save out plots
    mainFile = ['main_alpha_layer', num2str(layerIdx), '_stim1.png'];
    subFile = ['sub_alpha_layer', num2str(layerIdx), '_stim1.png'];
    saveas(mainFig, mainFile);
    saveas(subFig, subFile);
    
    % close figs
    close(mainFig);
    close(subFig);
                                          
end

                                  
                                  
                                  
                                  
                                  
                                  
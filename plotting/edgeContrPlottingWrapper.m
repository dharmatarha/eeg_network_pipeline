% % load necessary data: edge contributions and full raw connectivity
% e = load('~/eeg_network_pipeline/dev/p/edgeContributionsPlottingVars.mat');
% % only need raw connectivity and edge contributions as cohenD scores
% cohens = e.cohens;
% cohensM = e.cohensM;
% realConn = e.realConn;

% group level edge contributions and raw connectivity,
% alpha,
% iplv

e = load('/media/adamb/bonczData/hyperscan/newSurrEdgeEstimates/RS_vs_edgeContr_alpha_iplv.mat');
d = load('/media/adamb/bonczData/hyperscan/newSurrEdgeEstimates/alpha/group_surrResults_alpha_iplv.mat');

cohens = e.edgeDs; cohens(isnan(cohens)) = 0;
% cohens = mean(e.subEdgeDs,1)'; cohens(isnan(cohens)) = 0;
realConn = squeeze(mean(mean(d.meanConn, 2), 1));
groupRS = e.groupRS;
% matrix version of cohen d values
nodeNo = size(realConn, 1);
tmp = nan(nodeNo);
tmp(triu(true(nodeNo), 1)) = cohens;
cohensM = tmp;
realConnMean = realConn;

% % check if the matrix of contributions is arranged from the vector
% % correctly
% idx = tril(true(size(realConn, 1)), -1);
% c = nan(size(realConn, 1));
% c(idx) = cohens; c = c';
% if ~isequal(triu(c, 1), triu(cohensM, 1))
%     error('Bad arrangement of cohen D values from vector to matrix!');
% end

% get labels and colors
if ~exist('roiLabels', 'var')
    tmp = which('roiNamesInOrder.mat', '-ALL');
    if size(tmp, 1) ~= 1
        error('Found either zero or multiple versions of "roiNamesInOrder.mat"! See the help on input arg "roiLabels"!');
    else
        tmp = load('roiNamesInOrder.mat');
        roiLabels = tmp.roisShort;
    end
end
if ~exist('colorTriplets', 'var')
    tmp = which('colorTriplets.mat', '-ALL');
    if size(tmp, 1) ~= 1
        error('Found either zero or multiple versions of "colorTriplets.mat"! See the help on input arg "colorTriplets"!');
    else
        tmp = load('colorTriplets.mat');
        colorTriplets = tmp.colorTriplets24;
    end
end
% further checks
if length(roiLabels) ~= size(realConn, 1)
    error('Length of ROI / node labels cell array does not match the number of nodes in "realConn"!');
end
if ~iscolumn(roiLabels)
    roiLabels = roiLabels';
end

% define trimmingThr and group2color
trimmingThr = [0, 0.114];
group2color = [1, 1];

% % get averaged connectivity matrix
% realConnMean = mean(mean(realConn, 3, 'omitnan'), 4, 'omitnan');

% make sure we have symmetric matrices with zeros at diagonals before reordering rows/columns
realConnMean = triu(realConnMean, 1) + triu(realConnMean, 1)';
cohensM = triu(cohensM, 1) + triu(cohensM, 1)';

% reorder connectivity and edge grouping matrices + get new, correct labels
[equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(roiLabels);
if ~equalFlag && matchingSetsFlag
    % rearrange connectivity matrix based on new ROI label order
    [realConnReord, old2new1] = matrixReorder(realConnMean, roiLabels, roiLabelsPlotting);
    % rearrange edge grouping matrix based on new ROI label order
    [cohensMReord, old2new2] = matrixReorder(cohensM, roiLabels, roiLabelsPlotting);    
    % check for equality of rearranging vectors - sanity check
    if ~isequal(old2new1, old2new2)
        error('matrix rearranges do not match...');
    end
else
    error('Unexpected mismatch or equality in terms of labels !!!');
end

% define edge grouping matrix from the cohens values
cSort = sort(cohens, 'descend');
topN = 30;  % top contributing edges to highlight
edgeMembership = double(cohensMReord>cSort(topN+1));  % binary, only one group to highlight

% define final connectivity matrix with the right name
connMatrix = realConnReord;
% define final label cell array with the correct name
labels = roiLabelsPlotting;

% only upper triangles
connMatrix(tril(true(size(connMatrix, 1)))) = nan;
edgeMembership(tril(true(size(edgeMembership, 1)))) = nan;

% clear vars not needed anymore 
clearvars e old2new1 old2new2 tmp c idx

% call the main plotting function
mainFig = circleGraphPlot_edges(connMatrix, edgeMembership, colorTriplets,...
                                group2color, trimmingThr, labels, 'draw');




%% Version for circleGraphPlot_edgeColorWeights


% define edge grouping matrix from the cohens values
cSort = sort(cohens, 'descend');
topN = 60;  % top contributing edges to highlight
edgeMembership = double(cohensMReord>cSort(topN+1));  % binary, only one group to highlight

% define final connectivity matrix with the right name
connMatrix = realConnReord;
% define final label cell array with the correct name
labels = roiLabelsPlotting;

% only upper triangles
connMatrix(tril(true(size(connMatrix, 1)))) = nan;
edgeMembership(tril(true(size(edgeMembership, 1)))) = nan;

% filter for edges only in topN
connMatrixTrimmed = connMatrix;
connMatrixTrimmed(edgeMembership==0) = 0;

% coloring is based on cohensM
cohensMReordTrimmed = cohensMReord;
cohensMReordTrimmed(edgeMembership==0) = nan;

mainFig = circleGraphPlot_edgeColorWeights(connMatrixTrimmed, cohensMReordTrimmed, 'autumn',...
                                trimmingThr, labels, 'draw');

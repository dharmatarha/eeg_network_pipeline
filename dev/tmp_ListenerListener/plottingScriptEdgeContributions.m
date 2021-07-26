%% Load data

% load pruned group-level data
connFile = '/media/adamb/bonczData/hyperscan/alpha_iplv_GroupPruned_connData.mat';
tmp = load(connFile);
connData = tmp.connDataPruned;

% load edge contributions
edgeF = '/media/adamb/bonczData/hyperscan/edgeContributions/alpha_iplv_GroupPruned_cmpFullEdge_corr.mat';
tmp = load(edgeF);
ps = tmp.permRes.pEstAll;
ds = tmp.permRes.cohendAll;


%% Adjust p-values, derive a mask

% adjust p values
ps(ps==0) = 0.0001;
ps(ps<0.1) = ps(ps<0.1)*2;  % two-tailed

% % threshold - FDR
% q = 0.01;
% [h, crit] = fdr(ps, q);
% mask = ps<=crit;

% threshold - based on effect size
% dsMin = 0.119;
dsMin = 0.16;
mask = ds>dsMin;

disp(['No. of edges after thresholding: ', num2str(sum(mask))]);


%% Network similarity heatmap after masking (corr)

% prepare vars for heatmap
dataLin = linearizeTrius(connData);
dataLin = reshape(dataLin, [1891, 160]);
% mask edges
% dataLin(~mask, :) = [];
% correlate 
connSim = corr(dataLin);
% mask diagonal
connSim(1:161:end) = NaN;

% heatmap
hmap = heatmap(connSim, 'ColorMap', jet);
% remove tick labels - release specific method accessing private properties:
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(hmap);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];
hmap.FontSize = 14;
set(gcf, 'color', 'w');


%% Prepare vars for circle plot coloring edges based on contributions

% group avg conn matrix
avgConn = squeeze(mean(mean(connData, 4), 3));

% get mask matrix
maskM = nan(62);
maskM(triu(true(62), 1)) = mask;
maskM(maskM==0) = NaN;

% get ds matrix, masked
dsM = nan(62);
dsM(triu(true(62), 1)) = ds;
dsM(isnan(maskM)) = NaN;

% apply mask to conn matrix
connMasked = avgConn;
connMasked(isnan(maskM)) = 0;

% get ROIs
roiF = '/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat';
tmp = load(roiF);
labels = tmp.roisShort;
labels = labels';
% match with predefined sets
[equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(labels);
disp(['Matching rois? ', num2str(matchingSetsFlag)]);
if matchingSetsFlag
    % symmetrize matrices
    connMasked = triu(connMasked, 1)+triu(connMasked, 1)';
    dsM = triu(dsM, 1)+triu(dsM, 1)';
    dsM(1:63:end) = NaN;
    % reorder connectivity matrix
    [connMatrix, old2new] = matrixReorder(connMasked, labels, roiLabelsPlotting);
    % reorder ds matrix
    dsMatrix = dsM(old2new, old2new);
end

% select a colormap 
colorMap = 'autumn';

edgeColorWeights = dsMatrix;
trimmingThr = 0;
figTitle = 'Alpha, iplv';
labels = roiLabelsPlotting;

% call plotting function
mainFig = circleGraphPlot_edgeColorWeights(connMatrix, edgeColorWeights, colorMap, labels, trimmingThr, figTitle);


% %% Stats on connSim - from cmpFullConn
% 
% permNo = 10000;
% permStat = 'mean';
% epochNo = 40;
% condNo = 4;
% 
% % first mirror the upper triangular part of connSim matrix to the lower
% % half
% connSim = triu(connSim, 1) + triu(connSim, 1)'; 
% 
% % preallocate results struct
% permRes = struct;
% permRes.withinCondMean = nan;
% permRes.withinCondSD = nan;
% permRes.withinCondMedian = nan;
% permRes.acrossCondMean = nan;
% permRes.acrossCondSD = nan;
% permRes.acrossCondMedian = nan;
% permRes.realDiff = nan;
% permRes.permDiff = nan(permNo, 1);
% permRes.pEst = nan;
% 
% % preallocate vars holding within- and across condition/stimulus similarities
% withinCondSim = nan((epochNo^2-epochNo)/2, condNo);
% acrossCondSim = nan(epochNo^2*(condNo-1), condNo);
% 
% for condIdx = 1:condNo
%     
%     % within-cond epoch pairing similarities for given condition
%     tmp = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, (condIdx-1)*epochNo+1:condIdx*epochNo); 
%     withinCondSim(:, condIdx) = tmp(triu(true(epochNo), 1));
%     
%     % across-condition pairings for epochs in given condition
%     tmpParts = cell(2, 1);
%     if condIdx ~= 1
%         tmpParts{1} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, 1:(condIdx-1)*epochNo);
%     end
%     if condIdx ~= condNo
%         tmpParts{2} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, condIdx*epochNo+1:end);
%     end
%     % concatenate the parts
%     if isempty(tmpParts{1})
%         tmpAll = tmpParts{2};
%     elseif isempty(tmpParts{2})
%         tmpAll = tmpParts{1};
%     else
%         tmpAll = [tmpParts{1}, tmpParts{2}];
%     end
%     % linearize into a vector
%     acrossCondSim(:, condIdx) = reshape(tmpAll', [epochNo^2*(condNo-1), 1]);
%     
%     % store descriptives in result struct
%     permRes(condIdx).withinCondMean = mean(withinCondSim(:, condIdx));
%     permRes(condIdx).withinCondSD = std(withinCondSim(:, condIdx));
%     permRes(condIdx).withinCondMedian = median(withinCondSim(:, condIdx));
%     permRes(condIdx).acrossCondMean = mean(acrossCondSim(:, condIdx));
%     permRes(condIdx).acrossCondSD = std(acrossCondSim(:, condIdx)); 
%     permRes(condIdx).acrossCondMedian = median(acrossCondSim(:, condIdx));
%     
%     % permutation test across the two groups of connectivity similarity
%     % values
%     [permRes(condIdx).pEst,... 
%      permRes(condIdx).realDiff,... 
%      permRes(condIdx).permDiff,... 
%      permRes(condIdx).cohenD] = permTest(withinCondSim(:, condIdx),... 
%                                          acrossCondSim(:, condIdx),... 
%                                          permNo,... 
%                                          permStat,...
%                                          'silent');
%     
% end
% % user message
% disp([char(10), 'Compared similarity within- versus across-condition epoch-pairings']);
% 
% % Comparison between within- and across-condition pairings across all conditions / stimuli
% 
% % store descriptives in result struct
% permRes(condNo+1).withinCondMean = mean(withinCondSim(:));
% permRes(condNo+1).withinCondSD = std(withinCondSim(:));
% permRes(condNo+1).withinCondMedian = median(withinCondSim(:));
% permRes(condNo+1).acrossCondMean = mean(acrossCondSim(:));
% permRes(condNo+1).acrossCondSD = std(acrossCondSim(:)); 
% permRes(condNo+1).acrossCondMedian = median(acrossCondSim(:));
% 
% % permutation test on vectorized withinCondSim and acrossCondSim matrices
% [permRes(condNo+1).pEst,... 
%  permRes(condNo+1).realDiff,... 
%  permRes(condNo+1).permDiff,... 
%  permRes(condNo+1).cohenD] = permTest(withinCondSim(:),... 
%                                       acrossCondSim(:),... 
%                                       permNo,... 
%                                       permStat,...
%                                       'silent');
% % user message
% disp([char(10), 'Compared within- vs across-cond similarity across all stimuli']);    
%     
% 
% % Compare similarity values between different within-condition pairings
% % e.g. between epoch-pairings in within cond 1 and epoch-pairings within
% % cond 2
% 
% % preallocate results struct
% withinCondPermRes = struct;
% % loops through conditions
% counter = 0;
% for condOne = 1:condNo
%     for condTwo = 1:condNo
%         % only calculate if the two conditions are not the same
%         if condTwo > condOne
%             
%             % adjust counter
%             counter = counter + 1;
%             
%             % permutation test across the two groups of connectivity similarity
%             % values
%             [withinCondPermRes(counter).pEst,... 
%                 withinCondPermRes(counter).realDiff,... 
%                 withinCondPermRes(counter).permDiff,...
%                 withinCondPermRes(counter).cohenD] = permTest(withinCondSim(:, condOne),... 
%                                                               withinCondSim(:, condTwo),... 
%                                                               permNo,... 
%                                                               permStat,...
%                                                               'silent');
%             
%         end  % if 
%     end  % for condTwo
% end  % for condOne
% % user message
% disp([char(10), 'Compared similarity across conditions for within-condition epoch-pairings']);
% disp('Done with everything!');



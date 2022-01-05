%% Script for Figure 2A
% The script generates the group-averaged connectivity matrix, 
% reordered for easy-to-read labeling.
%
% 

%% Base params

method = 'ciplv'; 
freq = 'delta';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';

% option for jackknife-based thresholding of group mean data
jkknife = false; 
jkknife_threshold = 0.05;  % probability threshold that the sample mean is zero ( = level for confidence interval)

if ~strcmp(freq, 'delta')
    subNo = 26;
else
    subNo = 25;
end

condNo = 4;
epochNo = 40;
roiNo = 62;


%% Load data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['group_surrResultsv2_', freq, '_', method, '.mat']);
res = load(fileP);


%% Delete-one jackknife resampling if it was requested

if jkknife
    
    % preallocate jkknifed samples var
    jkSamples = nan([subNo, condNo, epochNo, roiNo, roiNo]);
    
    % load subject-level data
    switch thr
        case 'unthr'
            subConn = res.realConn;
        case 'thrSub'
            % get the data from the subject-level thresholded matrices, only positive diffs
            subConn = res.maskedConnPos;   
        case 'thrGroup'
            error('This is not really implemented, probably does not make sense either');
    end
    
%     % go through subjects, create delete-one means
%     for subIdx = 1:subNo
%         tmp = subConn; 
%         tmp(subIdx, :, :, :, :) = [];
%         jkSamples(subIdx, :, :, :, :) = mean(tmp, 1);
%     end
%     
%     
%     % test the samples for each edge
%     jk_hMap = nan(roiNo);  % map for edges passing / failing z-test
%     jk_pMap = jk_hMap;  % map for p values from z-tests
% 
%     for edge1 = 1:roiNo
%         for edge2 = edge1+1:roiNo
%             tmpSample = jkSamples(:, :, :, edge1, edge2);
%             tmpSample = tmpSample(:);  % cast into a vector
%             zeroP = normcdf(0, mean(tmpSample), std(tmpSample));
%             jk_hMap(edge1, edge2) = zeroP < jkknife_threshold;
%             jk_pMap(edge1, edge2) = zeroP;
%         end
%     end
% 
%     % get vector formats from maps
%     jk_h = jk_hMap(triu(true(roiNo), 1));
%     jk_p = jk_pMap(triu(true(roiNo), 1));


    % go through subjects, create delete-one means
    for subIdx = 1:subNo
        tmp = subConn; 
        tmp(subIdx, :, :, :, :) = [];
        jkSamples(subIdx, :, :, :, :) = mean(tmp, 1);
    end
    
    
    % test the samples for each edge
    jk_hMap = nan(roiNo);  % map for edges passing / failing z-test
    jk_pMap = jk_hMap;  % map for p values from z-tests

    for edge1 = 1:roiNo
        for edge2 = edge1+1:roiNo
            tmpSample = subConn(:, :, :, edge1, edge2);
            tmpSample = tmpSample(:);  % cast into a vector
            zeroP = ttest(tmpSample, 0);
            jk_hMap(edge1, edge2) = zeroP < jkknife_threshold;
            jk_pMap(edge1, edge2) = zeroP;
        end
    end

    % get vector formats from maps
    jk_h = jk_hMap(triu(true(roiNo), 1));
    jk_p = jk_pMap(triu(true(roiNo), 1));


%     % test the samples for each edge
%     jk_hMap = nan(condNo, epochNo, roiNo, roiNo);  % map for edges passing / failing z-test
%     jk_pMap = jk_hMap;  % map for p values from z-tests
% 
%     for condIdx = 1:condNo
%         for epochIdx = 1:epochNo
%             for edge1 = 1:roiNo
%                 for edge2 = edge1+1:roiNo
%                     tmpSample = subConn(:, condIdx, epochIdx, edge1, edge2);
%                     tmpSample = tmpSample(:);  % cast into a vector
%                     [~, zeroP] = ttest(tmpSample, 0);
%                     jk_hMap(condIdx, epochIdx, edge1, edge2) = zeroP < jkknife_threshold;
%                     jk_pMap(condIdx, epochIdx, edge1, edge2) = zeroP;
%                 end
%             end
%         end
%     end
% 
%     % get vector formats from maps
%     jk_h = nan(condNo, epochNo, roiNo*(roiNo-1)/2);
%     jk_p = jk_h;
%     for condIdx = 1:condNo
%         for epochIdx = 1:epochNo
%             tmp_hMap = squeeze(jk_hMap(condIdx, epochIdx, :, :));
%             jk_h(condIdx, epochIdx, :) = tmp_hMap(triu(true(roiNo), 1));
%             tmp_pMap = squeeze(jk_pMap(condIdx, epochIdx, :, :));
%             jk_p(condIdx, epochIdx, :) = tmp_pMap(triu(true(roiNo), 1));
%         end
%     end
    
end
            
    

%% Get mean connectivity, depending on "thr"
% Only versions where the positively larger connections are considered

switch thr
    case 'unthr'
        connData = res.meanConn;
    case 'thrSub'
        % get the data from the averaged, subject-level thresholded matrices, only positive diffs
        connData = res.meanMaskedConnPos;   
    case 'thrGroup'
        connData = res.groupMaskedConnPos;
end  % switch
        
% get averaged conn matrix
meanConn = squeeze(mean(mean(connData, 1), 2));


%% Load ROI label names

labels = load('/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;  % short names


%% Get reordered heatmap

[connMatrixRo, labelsRo, h, hs] = heatmapReordered(meanConn, labels);


%% Save out, quit

saveas(gcf, fullfile(baseDir, ['meanConn_', thr, '_', freq, '_', method, '.svg']));
saveas(gcf, fullfile(baseDir, ['meanConn_', thr, '_', freq, '_', method, '.png']));



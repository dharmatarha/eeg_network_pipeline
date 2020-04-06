% quick and dirty play with connectivity matrices


%% basic variables to determine what we work with

freq = 'alpha';
dirBase = ['/home/adamb/MATLAB/speaker_listener_data/', freq];
listDataF = [dirBase, '/', freq, '_angleConnectivity.mat'];
speakDataF = [dirBase, '/speaker_', freq, '_angleConnectivity.mat'];
 
% savefile for interim / final results
saveF = [dirBase, '/connectivityResults_', freq, '.mat'];
    
trialNo = 162;
%% load data

load(listDataF);
listPLI = pli; listPLV = plv;
load(speakDataF);
speakPLI = pli; speakPLV = plv;
clearvars pli plv;

% get sizes
listNo = size(listPLI, 1);
trialNo = size(listPLI, 2);
roiNo = size(listPLI, 3);


%% Linearize data

% preallocate
listPLIlin = nan(listNo, trialNo, roiNo*roiNo);
listPLVlin = nan(listNo, trialNo, roiNo*roiNo);
speakPLIlin = nan(trialNo, roiNo*roiNo);
speakPLVlin = nan(trialNo, roiNo*roiNo);

% listener
for l = 1:listNo
    for t = 1:trialNo
        listPLIlin(l, t, :) = reshape(listPLI(l, t, :, :), [roiNo*roiNo, 1]);
        listPLVlin(l, t, :) = reshape(listPLV(l, t, :, :), [roiNo*roiNo, 1]);
    end
end

% speaker
for t = 1:trialNo
    speakPLIlin(t, :) = reshape(speakPLI(t, :, :), [roiNo*roiNo, 1]);
    speakPLVlin(t, :) = reshape(speakPLV(t, :, :), [roiNo*roiNo, 1]);
end

% get rid of NaN values
% sanity check first
tmpList = squeeze(isnan(listPLIlin(1, 1, :)));
tmpSpeak = isnan(speakPLIlin(1, :));
if ~iscolumn(tmpSpeak)
    tmpSpeak = tmpSpeak';
end
if ~isequal(isnan(tmpList), isnan(tmpSpeak))
    error('NaN value locations do not match across linearized speaker and listener connectivity matrices');
end
% indices of NaN values
idx = find(tmpList);
% delete the slices
listPLIlin(:, :, idx) = [];
listPLVlin(:, :, idx) = [];
speakPLIlin(:, idx) = [];
speakPLVlin(:, idx) = [];

% get new size for connectivity vector
conNo = size(listPLIlin, 3);


%% Create thresholded version (only strngest connections)

% strongest x% survives, ratio for the threshold:
strRatio = 0.1;

% find thresholds
thrList = nan(listNo, trialNo, 2); % last dim is for PLI/PLV
thrSpeak = nan(trialNo, 2);
% listeners
for l = 1:listNo
    for t = 1:trialNo
        tmp = sort(squeeze(listPLIlin(l, t, :)));
        thrList(l, t, 1) = tmp(end-floor(conNo*strRatio));
        tmp = sort(squeeze(listPLVlin(l, t, :)));
        thrList(l, t, 2) = tmp(end-floor(conNo*strRatio)); 
    end
end
% speaker
for t = 1:trialNo
    tmp = sort(squeeze(speakPLIlin(t, :)));
    thrSpeak(t, 1) = tmp(end-floor(conNo*strRatio));
    tmp = sort(squeeze(speakPLVlin(t, :)));
    thrSpeak(t, 2) = tmp(end-floor(conNo*strRatio)); 
end
    
% use thresholds
% preallocate new variables
listPLIlinTH = listPLIlin;
listPLVlinTH = listPLVlin;
speakPLIlinTH = speakPLIlin;
speakPLVlinTH = speakPLVlin;
% listeners
for l = 1:listNo
    for t = 1:trialNo
        listPLIlinTH(l, t, listPLIlinTH(l, t, :)<thrList(l, t, 1)) = NaN;
        listPLVlinTH(l, t, listPLVlinTH(l, t, :)<thrList(l, t, 2)) = NaN;
    end
end
% speakers
for t = 1:trialNo
    speakPLIlinTH(t, speakPLIlinTH(t, :)<thrSpeak(t, 1)) = NaN;
    speakPLVlinTH(t, speakPLVlinTH(t, :)<thrSpeak(t, 2)) = NaN;
end


% %% Interim save
% 
% save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
%     'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
%     'thrList', 'thrSpeak');


%% Compare connectivity across listeners vs within-listeners

% calculate within-listener similarity
withinListPLI = nan(listNo, trialNo, trialNo);
withinListPLV = withinListPLI;

for l = 1:listNo    
    counter = 0;
    trialPairings = nan(trialNo*(trialNo-1)/2, 2);
    for t1 = 1:trialNo
        for t2 = 1:trialNo
            if t1~=t2 && ~ismember([t1, t2], trialPairings, 'rows') && ~ismember([t2, t1], trialPairings, 'rows')
                counter = counter+1;
                trialPairings(counter, :) = [t1, t2];
                withinListPLI(l, t1, t2) = corr(squeeze(listPLIlin(l, t1, :)), squeeze(listPLIlin(l, t2, :)));
                withinListPLV(l, t1, t2) = corr(squeeze(listPLVlin(l, t1, :)), squeeze(listPLVlin(l, t2, :)));
            end
        end
    end
    disp(['Done with listener ', num2str(l)]);
end


%% Interim save

save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
    'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
    'thrList', 'thrSpeak',...
    'withinListPLI', 'withinListPLV');


%% Across listeners
acrossListPLI = nan(listNo, listNo, trialNo);
acrossListPLV = acrossListPLI;

for l1 = 1:listNo
    for l2 = 1:listNo
        if l1~=l2 
            
            for t = 1:trialNo
                acrossListPLI(l1, l2, t) = corr(squeeze(listPLIlin(l1, t, :)), squeeze(listPLIlin(l2, t, :)));
                acrossListPLV(l1, l2, t) = corr(squeeze(listPLVlin(l1, t, :)), squeeze(listPLVlin(l2, t, :)));
            end
            
        end
        
    end
    
    disp(['Done with novel listener pairings for ', num2str(l1)]);
end


% %% Interim save
% 
% save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
%     'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
%     'thrList', 'thrSpeak',...
%     'withinListPLI', 'withinListPLV',...
%     'acrossListPLI', 'acrossListPLV');


%% Compare within- to across-

% within-listener
tmp = reshape(withinListPLI, listNo, trialNo^2);
withinListPLImeans = mean(tmp, 2, 'omitnan');
tmp = reshape(withinListPLV, listNo, trialNo^2);
withinListPLVmeans = mean(tmp, 2, 'omitnan');

% across-listener
tmp = nan(listNo, listNo*trialNo);
for l = 1:listNo
    tmp(l, :) = reshape(acrossListPLI(l, :, :), 1, listNo*trialNo);
end
acrossListPLImeans = mean(tmp, 2, 'omitnan');
tmp = nan(listNo, listNo*trialNo);
for l = 1:listNo
    tmp(l, :) = reshape(acrossListPLV(l, :, :), 1, listNo*trialNo);
end
acrossListPLVmeans = mean(tmp, 2, 'omitnan');


% %% Interim save
% 
% save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
%     'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
%     'thrList', 'thrSpeak',...
%     'withinListPLI', 'withinListPLV',...
%     'acrossListPLI', 'acrossListPLV',...
%     'withinListPLImeans', 'withinListPLVmeans',...
%     'acrossListPLImeans', 'acrossListPLVmeans');


%% First interesting plots

% mean of the within-listener connectivity matrices
meanWithinListPLI = squeeze(mean(withinListPLI, 1, 'omitnan'));
meanWithinListPLV = squeeze(mean(withinListPLV, 1, 'omitnan'));
% pli plot
heatmap(meanWithinListPLI);
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, [dirBase, '/acrossTrialsConnectivity_List_PLI.png']);
close;
% plv plot
heatmap(meanWithinListPLV);
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, [dirBase, '/acrossTrialsConnectivity_List_PLV.png']);
close;

% mean across-listener connectivity matrices - a peek at the variability
meanAcrossListPLI = squeeze(mean(acrossListPLI, 3, 'omitnan'));
meanAcrossListPLV = squeeze(mean(acrossListPLV, 3, 'omitnan'));
% pli plot
heatmap(meanAcrossListPLI);
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, [dirBase, '/acrossListenersConnectivity_PLI.png']);
close;
% plv plot
heatmap(meanAcrossListPLV);
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, [dirBase, '/acrossListenersConnectivity_PLV.png']);
close;



%% Connectivity similarity across speaker and listeners

% get averaged listener connectivities
meanListPLIlin = squeeze(mean(listPLIlin, 1, 'omitnan'));
meanListPLVlin = squeeze(mean(listPLVlin, 1, 'omitnan'));
meanListPLIlinTH = squeeze(mean(listPLIlinTH, 1, 'omitnan'));
meanListPLVlinTH = squeeze(mean(listPLVlinTH, 1, 'omitnan'));

% preallocate result matrices
speakListPLI = nan(trialNo, trialNo);
speakListPLV = nan(trialNo, trialNo);
speakListPLI_TH = nan(trialNo, trialNo);
speakListPLV_TH = nan(trialNo, trialNo);

% double loop through trials, get all correlations across connectivity
% matrices
for t1 = 1:trialNo
    for t2 = 1:trialNo
            
            speakListPLI(t1, t2) = corr(meanListPLIlin(t1, :)', speakPLIlin(t2, :)');
            speakListPLV(t1, t2) = corr(meanListPLVlin(t1, :)', speakPLVlin(t2, :)');
            speakListPLI_TH(t1, t2) = corr(meanListPLIlinTH(t1, :)', speakPLIlinTH(t2, :)', 'rows', 'complete');
            speakListPLV_TH(t1, t2) = corr(meanListPLVlinTH(t1, :)', speakPLVlinTH(t2, :)', 'rows', 'complete');            
        
    end
end



%% Interim save

save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
    'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
    'thrList', 'thrSpeak',...
    'withinListPLI', 'withinListPLV',...
    'acrossListPLI', 'acrossListPLV',...
    'withinListPLImeans', 'withinListPLVmeans',...
    'acrossListPLImeans', 'acrossListPLVmeans',...
    'meanWithinListPLI', 'meanWithinListPLV',...
    'meanAcrossListPLI', 'meanAcrossListPLV',...
    'speakListPLI', 'speakListPLV',...
    'speakListPLI_TH', 'speakListPLV_TH');


%% Speaker-listener comparisons: same-trial similarity vs corss-trial similarities

% concatenate all types of connectivity so we can easily loop through them
tmp = cat(3, speakListPLI, speakListPLV, speakListPLI_TH, speakListPLV_TH);

% preallocate results matrices
speakListPs = zeros(size(tmp, 3), 1);
speakListDiffs = zeros(size(tmp, 3), 1);

for repr = 1:size(tmp, 3)
    
    % prepare data for comparison: similarity from matching vs.
    % non.matching trials
    matchingTrials = diag(squeeze(tmp(:, :, repr)));
    nonmatchingTrials = squeeze(tmp(:, :, repr));
    nonmatchingTrials(1:size(nonmatchingTrials, 1)+1:end) = NaN;  % diagonal is filled with NaN values
    nonmatchingTrials = reshape(nonmatchingTrials, 1, size(nonmatchingTrials, 1)^2);
    nonmatchingTrials(isnan(nonmatchingTrials)) = [];
    
    % compare
    [speakListPs(repr), speakListDiffs(repr)] = permTest(matchingTrials, nonmatchingTrials, 10000);
    disp([char(10), 'Difference between matching (same trial) ',...
        'and non-matching (different trial) connectivity similarities: ',...
        num2str(speakListDiffs(repr)), char(10), 'Corresponding probability: ',...
        num2str(speakListPs(repr))]);
    % plot
    histogram(matchingTrials); 
    hold on; 
    histogram(nonmatchingTrials, 'DisplayStyle', 'stairs'); 
    legend('matching trials', 'non-matching trials'); 
    annotation('textbox', [0.2, 0.5, 0.3, 0.3],... 
        'String', ['p = ', num2str(speakListPs(repr))],... 
        'FitBoxToText','on');
    hold off;
    % save figure as png
    set(gcf, 'Position', get(0, 'Screensize'))
    saveas(gcf, [dirBase, '/speakerListener_connectivitySimilarity_', num2str(repr), '.png']);
    close;
    
end


%% Interim save

save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
    'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
    'thrList', 'thrSpeak',...
    'withinListPLI', 'withinListPLV',...
    'acrossListPLI', 'acrossListPLV',...
    'withinListPLImeans', 'withinListPLVmeans',...
    'acrossListPLImeans', 'acrossListPLVmeans',...
    'meanWithinListPLI', 'meanWithinListPLV',...
    'meanAcrossListPLI', 'meanAcrossListPLV',...
    'speakListPLI', 'speakListPLV',...
    'speakListPLI_TH', 'speakListPLV_TH',...
    'speakListPs', 'speakListDiffs');



%% Speaker-listener comparisons with +1 lag: +1 trial similarity vs corss-trial similarities

% concatenate all types of connectivity so we can easily loop through them
tmp = cat(3, speakListPLI, speakListPLV, speakListPLI_TH, speakListPLV_TH);

% preallocate results matrices
speakListPs_plus1 = zeros(size(tmp, 3), 1);
speakListDiffs_plus1 = zeros(size(tmp, 3), 1);

for repr = 1:size(tmp, 3)
    
    % prepare data for comparison: similarity from matching vs.
    % non.matching trials
    matchingTrials = squeeze(tmp(:, :, repr));
    matchingTrials = matchingTrials(2:size(matchingTrials, 1)+1:end);
    nonmatchingTrials = squeeze(tmp(:, :, repr));
    nonmatchingTrials(1:size(nonmatchingTrials, 1)+1:end) = NaN;  % diagonal is filled with NaN values
    nonmatchingTrials(2:size(nonmatchingTrials, 1)+1:end) = NaN;  % +1 lag is filled with NaN values
    nonmatchingTrials = reshape(nonmatchingTrials, 1, size(nonmatchingTrials, 1)^2);
    nonmatchingTrials(isnan(nonmatchingTrials)) = [];
    
    % compare
    [speakListPs_plus1(repr), speakListDiffs_plus1(repr)] = permTest(matchingTrials, nonmatchingTrials, 10000);
    disp([char(10), 'Difference between matching (same trial) ',...
        'and non-matching (different trial) connectivity similarities: ',...
        num2str(speakListDiffs_plus1(repr)), char(10), 'Corresponding probability: ',...
        num2str(speakListPs_plus1(repr))]);
    % plot
    histogram(matchingTrials); 
    hold on; 
    histogram(nonmatchingTrials, 'DisplayStyle', 'stairs'); 
    legend('matching trials', 'non-matching trials'); 
    annotation('textbox', [0.2, 0.5, 0.3, 0.3],... 
        'String', ['p = ', num2str(speakListPs_plus1(repr))],... 
        'FitBoxToText','on');
    hold off;
    % save figure as png
    set(gcf, 'Position', get(0, 'Screensize'))
    saveas(gcf, [dirBase, '/speakerListener_connectivitySimilarity_', num2str(repr), '_plus1.png']);
    close;
    
end


%% Interim save

save(saveF, 'listPLIlin', 'listPLVlin', 'listPLIlinTH', 'listPLVlinTH',...
    'speakPLIlin', 'speakPLVlin', 'speakPLIlinTH', 'speakPLVlinTH',...
    'thrList', 'thrSpeak',...
    'withinListPLI', 'withinListPLV',...
    'acrossListPLI', 'acrossListPLV',...
    'withinListPLImeans', 'withinListPLVmeans',...
    'acrossListPLImeans', 'acrossListPLVmeans',...
    'meanWithinListPLI', 'meanWithinListPLV',...
    'meanAcrossListPLI', 'meanAcrossListPLV',...
    'speakListPLI', 'speakListPLV',...
    'speakListPLI_TH', 'speakListPLV_TH',...
    'speakListPs', 'speakListDiffs',...
    'speakListPs_plus1', 'speakListDiffs_plus1');



%% Speaker-speaker connectivity similarity

speakWithinPLI = nan(trialNo, trialNo);
speakWithinPLV = nan(trialNo, trialNo);
speakWithinPLI_TH = nan(trialNo, trialNo);
speakWithinPLV_TH = nan(trialNo, trialNo);

for t1 = 1:trialNo
    for t2 = 1:trialNo
        if t1~=t2
            speakWithinPLI(t1, t2) = corr(speakPLIlin(t1, :)', speakPLIlin(t2,:)');
            speakWithinPLV(t1, t2) = corr(speakPLVlin(t1, :)', speakPLVlin(t2,:)');
            speakWithinPLI_TH(t1, t2) = corr(speakPLIlinTH(t1, :)', speakPLIlinTH(t2,:)', 'rows', 'complete');
            speakWithinPLV_TH(t1, t2) = corr(speakPLVlinTH(t1, :)', speakPLVlinTH(t2,:)', 'rows', 'complete');
        end
    end
end
    
    
    
%% CLose / Far
winL = 5;

% far
speakList_farCorr = zeros(10,4);

counter = 0;
for tS = 1:trialNo
    
    wintS = tS-winL:1:tS+winL;
    ts = 1:trialNo;
    ts(ismember(ts, wintS)) = [];
    tsL = max(size(ts));
    
    for tL = 1:tsL
        counter = counter+1;
        tCurrentL = ts(tL);
        
        speakList_farCorr(counter, 1) = corr(speakPLIlin(tS,:)', meanListPLIlin(tCurrentL,:)');
        speakList_farCorr(counter, 2) = corr(speakPLVlin(tS,:)', meanListPLVlin(tCurrentL,:)');
        speakList_farCorr(counter, 3) = corr(speakPLIlinTH(tS,:)', meanListPLIlinTH(tCurrentL,:)', 'rows', 'complete');
        speakList_farCorr(counter, 4) = corr(speakPLVlinTH(tS,:)', meanListPLVlinTH(tCurrentL,:)', 'rows', 'complete');
        
    end
end
    

% close
speakList_closeCorr = zeros(10,4);

counter = 0;
for tS = 1:trialNo
    
    wintS = tS-winL:1:tS+winL;
    ts = 1:trialNo;
    ts(~ismember(ts, wintS)) = [];
    tsL = max(size(ts));
    
    for tL = 1:tsL
        counter = counter+1;
        tCurrentL = ts(tL);
        
        speakList_closeCorr(counter, 1) = corr(speakPLIlin(tS,:)', meanListPLIlin(tCurrentL,:)');
        speakList_closeCorr(counter, 2) = corr(speakPLVlin(tS,:)', meanListPLVlin(tCurrentL,:)');
        speakList_closeCorr(counter, 3) = corr(speakPLIlinTH(tS,:)', meanListPLIlinTH(tCurrentL,:)', 'rows', 'complete');
        speakList_closeCorr(counter, 4) = corr(speakPLVlinTH(tS,:)', meanListPLVlinTH(tCurrentL,:)', 'rows', 'complete');
        
    end
end



%% other

% mean speaker connectivity matrix, based on PLV
meanSpeakConn = squeeze(mean(speakPLV, 1, 'omitnan'));
% thresholded version
targetRatio = 0.1;  
linTmp = reshape(meanSpeakConn, 1, size(meanSpeakConn,1)^2);
linTmp(isnan(linTmp)) = [];
sortedTmp = sort(linTmp);
thr = sortedTmp(floor(length(sortedTmp)*(1-targetRatio)));
meanSpeakConnThr = meanSpeakConn;
meanSpeakConnThr(meanSpeakConnThr<thr) = NaN;

% mean listener connectivity matrix, based on PLV
meanListConn = squeeze(mean(mean(listPLV, 1, 'omitnan'), 2, 'omitnan'));
% thresholded version
targetRatio = 0.1;  
linTmp = reshape(meanListConn, 1, size(meanListConn,1)^2);
linTmp(isnan(linTmp)) = [];
sortedTmp = sort(linTmp);
thr = sortedTmp(floor(length(sortedTmp)*(1-targetRatio)));
meanListConnThr = meanListConn;
meanListConnThr(meanListConnThr<thr) = NaN;

% shared connectivity
sharedConn = ~isnan(meanListConnThr) & ~isnan(meanSpeakConnThr);
% numeric version for viewing
sharedConnNum = double(sharedConn);
% averaged connectivity values
tmpList = meanListConnThr;
tmpList(~sharedConn) = NaN;
tmpSpeak = meanSpeakConnThr;
tmpSpeak(~sharedConn) = NaN;
avgSharedConn = (tmpList + tmpSpeak)/2;


%% look for maximum correlation in time between speaker and listener

% we go through all roi-roi connections, both in listener and speaker, and look for maximum
% correlation regions in time across speaker and listener
% based on PLV
meanListPLV = squeeze(mean(listPLV, 1, 'omitnan'));

% overkill version
corrRes = nan(roiNo, roiNo, roiNo, roiNo);
for speakRoi1 = 1:roiNo
    for speakRoi2 = 1:roiNo
        for listRoi1 = 1:roiNo
            for listRoi2 = 1:roiNo
                if speakRoi1~=speakRoi2 && listRoi1~=listRoi2
                
                    corrRes(speakRoi1, speakRoi2, listRoi1, listRoi2) = corr(squeeze(speakPLV(:, speakRoi1, speakRoi2)), squeeze(meanListPLV(:, listRoi1, listRoi2)));
                
                end
            end
        end
    end
end

% normal version
corrRes = nan(roiNo, roiNo);
for roi1 = 1:roiNo
    for roi2 = 1:roiNo
        corrRes(roi1, roi2) = corr(squeeze(meanListPLV(:, roi1, roi2)), squeeze(speakPLV(:, roi1, roi2)));
    end
end

% thresholding
targetRatio = 0.1;  
linTmp = reshape(corrRes, 1, size(corrRes,1)^2);
linTmp(isnan(linTmp)) = [];
sortedTmp = sort(linTmp);
thr = sortedTmp(floor(length(sortedTmp)*(1-targetRatio)));
corrResThr = corrRes;
corrResThr(corrResThr<thr) = NaN;

% logical version of thresholded map
corrResThrBoolean = corrResThr;
corrResThrBoolean(isnan(corrResThrBoolean)) = 0;
corrResThrBoolean = boolean(corrResThrBoolean);

% intersection of masks
speakerListenerCorrelated = corrResThrBoolean&sharedConn;


    




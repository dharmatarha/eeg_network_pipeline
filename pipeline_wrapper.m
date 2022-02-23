% Assume that both real and surrogate data are available, e.g. lv_s17_delta_plv.mat
% for real data and lv_s17_delta_surrEdgeEstReal_plv_iplv_ampCorr.mat for
% surrogate data, for each 202 tested person in each frequency band (delta, theta,
% alpha, beta, gamma).
%
% selectEpochsWrapper.m has to be run. This script calls selectEpochs() and
% sortSurrConn(). As a result, 75 non-overlapping epochs will be selected
% for each participant, and surrogate connectivity distribution parameters will
% be assigned to the selected epochs of the real connectivity data. The indices of
% the 75 selected epochs were randomly chosen for each participant, but the 
% indices are fixed for frequency bands. So random indices are selected for
% the first frequency band, and after that, these indices are used for all
% other frequency bands. This step makes result from different frequency
% bands comparable. Resulting files e.g. delta_plv_group.mat for selected
% real epochs and delta_plv_surrConn.mat for assigned surrogate connectivity
% distribution parameters. Moreover, in the surrogate file
% (e.g. delta_plv_surrConn.mat), there is thresholded data, resulting from
% both epoch-based thresholding, and also thresholding across epochs. In
% the further steps, thresholded data from epoch-based thresholding will not 
% be used because of the high pruning ration. Only thresholded data from
% thresholding across epochs will be used.
%
% stats/withinSubjectCorrelationNoWrap_thresholded() has to be run for all five
% frequency bands and all four connectivity metrics separately. This function
% calls connSimTest_subject_thresholding() wich calculates the
% within-subject reliability of thresholded data, in 1000 iterations for each
% participant. The resulting files e.g. delta_plv_simResCorr_thr.mat
% contains the 200 X 1000 within-subject reliability values for thresholded
% data. Using these files, violin plots can be created.
%
% stats/withinSubjectCorrelationWrapper_unthresholded.m has to be run. This
% script calls connSimTest_subject() for each connectivity metric and each
% frequency band. connSimTest_subject() calculates the within-subject 
% reliability of unthresholded data, in 1000 iterations for each participant. 
% The resulting files e.g. delta_plv_simResCorr_unthr.mat contains the
% 200 X 1000 within-subject reliability values for unthresholded data.
% Using these data, violin plots are also created by the script. The reason
% for using a wrapper for unthresholded data and using separate calls in
% each freq. and method for thresholded data is that thresholding is a time
% consuming procedure and a wrapper for thresholded data would run for more
% days.
%
% stats/groupSizeDependentCorrelationWrapper_thr.m has to be run. This
% script calls connSimTest_group_var_size_epochAveraged() and creates
% violin plots for each freqency band and connectivity metric, where
% between group correlation is calculated as a function of group size, for
% thresholded data (using the output of sortSurrConn()) with 1000 iterations
% for each group size, group sizes between 1 and 100.
%
% stats/groupSizeDependentCorrelationWrapper_unthr.m has to be run. This
% script calls connSimTest_group_var_size() and creates violin plots for 
% each freqency band and connectivity metric, where between group correlation
% is calculated as a function of group size, for unthresholded data with 1000 
% iterations for each group size, group sizes between 1 and 100.
% 
% stats/edgeContrToReliability_thr.m has to be run. This script calls 
% connSimTest_group_edgeContr_thr() and creates group-level average circle
% plots for each freqency band and connectivity metric, for thresholded
% data (using the output of sortSurrConn()). In the circle plot, edge width
% corresponds to connectivity strength and edge color corresponds to edge
% contribution to the reliability.
% 
% stats/edgeContrToReliability_unthr.m has to be run. This script calls 
% connSimTest_group_edgeContr() and creates group-level average circle
% plots for each freqency band and connectivity metric, for unthresholded
% data. In the circle plot, edge width corresponds to connectivity strength
% and edge color corresponds to edge contribution to the reliability.
% 
% epochDistanceSimilarityWrapper.m has to be run. This script calls 
% epochDistSim_subject() for each frequency band and each connectivity
% metric. As a result, the first 75 non-overlapping epochs are selected
% for each participant, and the similarity between all possible pairings of
% the selected epochs is calculated as a function of the distance between the
% paired epochs. This is calculated only for unthresholded data. Resulting
% files e.g. delta_plv_epochDistSim.mat. Results can be plotted using the 
% epochDistSimHelper.m script which plots the group-level average epoch
% similarity values for each epoch distance, each frequency band and each
% connectivity metric.


% run selectEpochsWrapper.m
% withinSubjectCorrelationNoWrap_thresholded('delta', 'plv')
% withinSubjectCorrelationNoWrap_thresholded('delta', 'iplv')
% withinSubjectCorrelationNoWrap_thresholded('delta', 'ampCorr')
% withinSubjectCorrelationNoWrap_thresholded('delta', 'orthAmpCorr')
% withinSubjectCorrelationNoWrap_thresholded('theta', 'plv')
% withinSubjectCorrelationNoWrap_thresholded('theta', 'iplv')
% withinSubjectCorrelationNoWrap_thresholded('theta', 'ampCorr')
% withinSubjectCorrelationNoWrap_thresholded('theta', 'orthAmpCorr')
% withinSubjectCorrelationNoWrap_thresholded('alpha', 'plv')
% withinSubjectCorrelationNoWrap_thresholded('alpha', 'iplv')
% withinSubjectCorrelationNoWrap_thresholded('alpha', 'ampCorr')
% withinSubjectCorrelationNoWrap_thresholded('alpha', 'orthAmpCorr')
% withinSubjectCorrelationNoWrap_thresholded('beta', 'plv')
% withinSubjectCorrelationNoWrap_thresholded('beta', 'iplv')
% withinSubjectCorrelationNoWrap_thresholded('beta', 'ampCorr')
% withinSubjectCorrelationNoWrap_thresholded('beta', 'orthAmpCorr')
% withinSubjectCorrelationNoWrap_thresholded('gamma', 'plv')
% withinSubjectCorrelationNoWrap_thresholded('gamma', 'iplv')
% withinSubjectCorrelationNoWrap_thresholded('gamma', 'ampCorr')
% withinSubjectCorrelationNoWrap_thresholded('gamma', 'orthAmpCorr')
% run withinSubjectCorrelationWrapper_unthresholded.m
% run groupSizeDependentCorrelationWrapper_thr.m
% run groupSizeDependentCorrelationWrapper_unthr.m
% run edgeContrToReliability_thr.m
% run epochDistanceSimilarityWrapper.m



% dirName = '../NAS502/EEG_resting_state/delta/';
% connMeasure = 'ampCorr';
% frequencyBand = 'delta';
% thresholding = 'thr';
% graphDistMetric = 'adjacencySpectral';
% surrType = 'edgesRandom';
% surrNo = 10;
% 
% mdsCorrWrapper(dirName, connMeasure, frequencyBand, thresholding, graphDistMetric, surrType, surrNo);
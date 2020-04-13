function prunedMask = getPrunedMask(prunedData, mode, thrs)
%% Calculates a mask of edges for a set of pruned connectivity matrices
%
% USAGE: prunedMask = getPrunedMask(prunedData, mode, thrs)
%
% A set of pruned connectivity matrices - even if from the same subject and
% task - is usually composed of partially overlapping elements. This 
% function derives a common mask based on the frequency of occurance of
% edges across the whole set.
% As our use cases are EEG datasets with size [ROIs, ROIs, epochs,
% conditions], inputs are either 4D or 5D arrays. In the latter case, the
% 5th dimension is for subjects. Masks can be either 2D (size [ROIs, ROIs])
% or 4D (size [ROIs, ROIs, epochs, conditions], only in the case of a 
% 5D input array). Arg "mode" controls which type of output is calculated.
% Edge selection in the mask is controlled by simple thresholds (arg
%

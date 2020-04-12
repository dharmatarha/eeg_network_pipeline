function [permRes, withinSim, acrossSim] = cmp2GroupMeanPruned(groupData, subIdx, varargin)
%% Compare individual connectivity to leave-one-out group mean connectivity
% Version for pruned (sparse) connectivity matrices.
%
% USAGE [permRes, withinSim, acrossSim] = cmp2GroupMeanPruned(groupData, subIdx, metric='corr', permNo=10000, permTest='mean')
%
% Function to compare the pruned connectivity matrix of a subject to that
% of the whole group except that subject (leave-one-out mean). 
% We evaluate the null hypothesis that the similarity of individual and
% group connectivity in matching epoch-pairings is the same as in 
% non-matching epoch-pairings.
% That is, similarity calculated from epochs e.g. 1-1, 2-2, ..., 
% 10-10, 11-11,... for the individual
% and the group, respectively, is the same as from epochs e.g. 1-30, 2-5,
% etc.
%
% We calculate both groups of similarity values (within-epoch vs
% across-epochs) and compare the two sets of values with a permutation test
%
% Mandatory inputs:
% groupData     - Numeric array, sized [ROIs, ROIs, epochs, conditions, subjects]. 
%           Group-level pruned connectivity matrix. Undirected
%           connectivity, so we work only with values in the upper
%           triangle of each connectivity matrix (connectivity matrix for 
%           a given epoch, condition and subject is defined by 
%           dimensions 1 and 2).
% subIdx        - Numeric value. Index of the subject whose data we compare 
%           with the rest of the group. Obviously needs to be one of
%           1:numberOfSubjects.
%
% Optional inputs:
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these distances is used: {'corr', 'eucl'}.  
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6.
% permStat      - String specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m
%
% Outputs:
% permRes       - Struct containing the output of permTest.m. Its fields
%               are pEst (estimated prob.), realDiff (the real difference
%               between within-epoch and across-epoch pairing similarities
%               in terms of the test stat), permDiff (the permuted
%               difference values) and CohenD (effect size, Cohen's d).
% withinSim     - Numeric vector of within-epoch pairing connectivity
%               similarities. Its size is
%               [numberOfEpochs*numberOfConditions, 1]
% acrossSim     - Numeric vector of across-epoch pairing connectivity
%               similarities. Its size is
%               [(numberOfEpochs^2*numberOfConditions^2-numberOfEpochs*numberOfConditions)/2, 1]
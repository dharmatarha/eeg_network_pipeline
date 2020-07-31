function [consensus, consensus_simm, pairwise_simm] = consensus_similarity(C)
%CONSENSUS_SIMILARITY     Construct a consensus (representative) partition
% using the iterative thresholding procedure
%
%   [consensus, consensus_simm, pairwise_simm] = CONSENSUS_SIMILARITY(C)
%   identifies a single representative partition from a set of C partitions
%   that is the most similar to the all others. Here, similarity is taken
%   to be the z-score of the Rand coefficient (see zrand.m)
%
%   NOTE: This code requires zrand.m to be on the MATLAB path
%
%   Inputs:     C,      pxn matrix of community assignments where p is the
%                       number of optimizations and n the number of nodes
%
%   Outputs:    consensus,      consensus partition
%               consensus_simm,	average similarity between consensus
%                               partition and all others
%               pairwise_simm,	pairwise similarity matrix
%   _______________________________________________
%   Marcelo G Mattar (08/21/2014) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downloaded from the Network Community Toolbox in 07/2020:
% http://commdetect.weebly.com/
%
% Reference to cite:
% Danielle S. Bassett, Mason A. Porter, Nicholas F. Wymbs, Scott T. Grafton, 
% Jean M. Carlson, Peter J. Mucha. Robust detection of dynamic community 
% structure in networks. Chaos, 2013, 23, 1. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npart = numel(C(:,1)); % number of partitions

% Initialize variables
pairwise_simm = zeros(npart,npart);

%% CALCULATE PAIRWISE SIMILARITIES
for i=1:npart
    for j=(i+1):npart
        pairwise_simm(i,j) = zrand(C(i,:),C(j,:));
    end
end
pairwise_simm = pairwise_simm + pairwise_simm';

% Average pairwise similarity
average_pairwise_simm = sum(pairwise_simm,2)/(npart-1);

%% EXTRACT PARTITION MOST SIMILAR TO THE OTHERS
[X,I] = max(average_pairwise_simm);
consensus = C(I,:);
consensus_simm = X;

return

function C=clustering_coef_wu(W)
%CLUSTERING_COEF_WU     Clustering coefficient
%
%   C = clustering_coef_wu(W);
%
%   The weighted clustering coefficient is the average "intensity"
%   (geometric mean) of all triangles associated with each node.
%
%   Input:      W,      weighted undirected connection matrix
%                       (all weights must be between 0 and 1)
%
%   Output:     C,      clustering coefficient vector
%
%   Note:   All weights must be between 0 and 1.
%           This may be achieved using the weight_conversion.m function,
%           W_nrm = weight_conversion(W, 'normalize');
%
%   Reference: Onnela et al. (2005) Phys Rev E 71:065103
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downloaded from the Brain Connectivity Toolbox in 10/2020:
% http://www.brain-connectivity-toolbox.net
%
% Reference to cite:
% Complex network measures of brain connectivity: Uses and interpretations.
% Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

K=sum(W~=0,2);            	
cyc3=diag((W.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
C=cyc3./(K.*(K-1));         %clustering coefficient

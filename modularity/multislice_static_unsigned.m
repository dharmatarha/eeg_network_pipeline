function [S Q lAlambda] = multislice_static_unsigned(A,gplus)

% INPUTS
% A is the (weighted) connectivity matrix
% it is assumsed that all values of the connectivity matrix are positive
% Gplus is the resolution parameter. If unsure, use default value of 1.
%
% OUTPUTS
% S is the partition (or community assignment of all nodes to communities)
% Q is the modularity of the (optimal) partition
% lAlambda is the effective fraction of antiferromagnetic edges (see Onnela
% et al. 2011 http://arxiv.org/pdf/1006.5731v1.pdf)
%
% This code uses the Louvain heuristic
%
% DB 2012

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


Aplus=A; Aplus(A<0)=0;
kplus=sum(Aplus)';
P = (kplus*kplus'/sum(kplus));
B=A-P.*gplus;
lAlambda = numel(find((A./P)<gplus));
%[S, Q] = greedysparse(B);
[S Q] = genlouvain(B, 10000, 0, 1, 'moverandw');
%[S2,Q2]=newmanklB(S,B); IF YOU WANT ITERATIVE IMPROVEMENT
Q=Q/(sum(kplus));




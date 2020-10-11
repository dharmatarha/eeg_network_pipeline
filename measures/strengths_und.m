function [str] = strengths_und(CIJ)
%STRENGTHS_UND        Strength
%
%   str = strengths_und(CIJ);
%
%   Node strength is the sum of weights of links connected to the node.
%
%   Input:      CIJ,    undirected weighted connection matrix
%
%   Output:     str,    node strength
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

% compute strengths
str = sum(CIJ);        % strength



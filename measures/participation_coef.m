function P=participation_coef(W,Ci,flag)
%PARTICIPATION_COEF     Participation coefficient
%
%   P = participation_coef(W,Ci);
%
%   Participation coefficient is a measure of diversity of intermodular
%   connections of individual nodes.
%
%   Inputs:     W,      binary/weighted, directed/undirected connection matrix
%               Ci,     community affiliation vector
%               flag,   0, undirected graph (default)
%                       1, directed graph: out-degree
%                       2, directed graph: in-degree
%
%   Output:     P,      participation coefficient
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
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

if ~exist('flag','var')
    flag=0;
end

switch flag
    case 0 % no action required
    case 1 % no action required
    case 2; W=W.';
end

n=length(W);                        %number of vertices
Ko=sum(W,2);                        %degree
Gc=(W~=0)*diag(Ci);                 %neighbor community affiliation
Kc2=zeros(n,1);                     %community-specific neighbors

for i=1:max(Ci)
   Kc2=Kc2+(sum(W.*(Gc==i),2).^2);
end

P=ones(n,1)-Kc2./(Ko.^2);
P(~Ko)=0;                           %P=0 if for nodes with no (out)neighbors

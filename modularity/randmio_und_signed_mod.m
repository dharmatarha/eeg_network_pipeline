function [R,eff] = randmio_und_signed_mod(W, ITER)
% RANDMIO_UND_SIGNED	Random graph with preserved signed degree distribution
%
%   R       = randmio_und_signed(W,ITER);
%   [R,eff] = randmio_und_signed(W,ITER);
%
%   This function randomizes an undirected network with positively and
%   negatively signed connections, while preserving the positively and
%   negatively signed degree distribution. The function does not preserve
%   the strength distribution in weighted networks.
%
%   Input:      W,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Downloaded from the Network Community Toolbox in 08/2020:
%   https://sites.google.com/site/bctnet/null
%
%   Reference:  Maslov and Sneppen (2002) Science 296:910
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin('randperm')==1
    warning('This function requires a recent (>2011) version of MATLAB.')
end

R     = double(W);              % sign function requires double input
n     = size(R,1);
ITER  = ITER*n*(n-1)/2;

% maximal number of rewiring attempts per 'iter'
maxAttempts = round(n/2);
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)    %while not rewired
        %select four distinct vertices
        nodes = randperm(n,4);
        a = nodes(1);
        b = nodes(2);
        c = nodes(3);
        d = nodes(4);
        
        r0_ab = R(a,b);
        r0_cd = R(c,d);
        r0_ad = R(a,d);
        r0_cb = R(c,b);
        
        %rewiring condition
        if (~isnan(r0_ab)) && (~isnan(r0_cd)) && (~isnan(r0_ad)) && (~isnan(r0_cb))
            
            R(a,d)=r0_ab;
            R(a,b)=r0_ad;
            R(c,b)=r0_cd;
            R(c,d)=r0_cb;
            
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations
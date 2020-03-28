function [prunedMatrix] = pruningFunction(pValuesMatrix, matrixToPrune, q)
%% Edge pruning based on FDR on p-value matrix
% 
% USAGE: prunedMatrix = pruningFunction(pValuesMatrix, matrixToPrune, q=0.05)
%
% Prunes the input matrix by first calculating FDR on the p-values and 
% then setting nonsignificant elements of the connectivity matrix to NaN.
%
% Mandatory inputs: 
% pValuesMatrix     - Matrix of p-values corresponding to the values in 
%                   matrixToPrune
% matrixToPrune     - Connectivity matrix where each value corresponds to 
%                   the value of an edge.
%
% Optional input:
% q                 - Q for FDR rate, numeric between 10^-9 - 1. 
%                   Defaults to 0.05. 
%
% Output:
% prunedMatrix     - Pruned connectivity matrix.
%


%% Input checks

% number of input args
if ~ismember(nargin, [2 3]) 
    error(['Function pruningFunction requires mandatory input args ',...
        '"pValuesMatrix" and "matrixToPrune" and optional arg "q"!']);
end
% set / check q
if nargin == 2
    q = 0.05;
else
    if ~isnumeric(q) || q>1 || q<10^-9 
        error('Optional input arg "q" is expected to be a number between 10e-9 - 1!');
    end
end
% matrix?
if ~ismatrix(pValuesMatrix) || ~ismatrix(matrixToPrune)
    error('Input args "pValuesMatrix" and "matrixToPrune" should be matrices!');
end
% check the size of input matrices
if ~isequal(size(pValuesMatrix), size(matrixToPrune))
    error('Input marices "pValuesMatrix" and "matrixToPrune" have different sizes!');
end
% sanity check: square matrices?
if ~isequal(size(pValuesMatrix, 1), size(pValuesMatrix, 1))
    error('Input args "pValuesMatrix" and "matrixToPrune" should be square matrices!')
end


%% FDR on p-values

% reshape and if there are NaN values, delete them
pValuesVect = reshape(pValuesMatrix, [size(pValuesMatrix, 1)^2, 1]);
pValuesVect(isnan(pValuesVect)) = [];

% FDR
[~, pCrit] = fdr(pValuesVect, q);


%% Set non-significant values in matrixToPrune to NaN

prunedMatrix = matrixToPrune;
prunedMatrix(pValuesMatrix>pCrit) = nan;


return
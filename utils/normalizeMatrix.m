function [normalizedMatrix] = normalizeMatrix(inputMatrix)
%% Matrix normalization
% 
% USAGE: normalizedMatrix = normalizeMatrix(inputMatrix)
%
% Matrix normalization is carried out such that the sum of all elements in the normalized matrix should be one.
%
% Mandatory input:
% inputMatrix       - The matrix to be normalized
%
% Output:
% normalizedMatrix  - Normalized matrix
%


%% Input checks

% Check if the input is a matrix
if ~ ismatrix(inputMatrix)
    error('Input shoud be a matrix');
end


%% Calculation

if any(any(isnan(inputMatrix)))
    warning('There were NaN values in the input matrix. We use ''omitnan'' flags, so the result will be correct.');
end
sumOfAllElements = sum(sum(inputMatrix, 'omitnan'), 'omitnan');
normalizedMatrix = inputMatrix / sumOfAllElements;



return
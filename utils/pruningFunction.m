function [prunedMatrix] = pruningFunction(pValuesMatrix, matrixToPrune)
%% Execute pruning
% 
% USAGE: prunedMatrix = pruningFunction(pValuesMatrix, matrixToPrune)
%
% Prunes the input matrix based on the p-values by setting nonsignificant
% matrix elements to NaN.
%
% Input: 
% pValuesMatrix     - Matrix of p-values corresponding to members in the matrix
%              to prune.
% matrixToPrune    - Matrix to prune.
%
% Output:
% prunedMatrix     - Pruned matrix.
%
% Notes: It is assumed that the input matrix "matrixToPrune" is symmetrical
% and pruning is executed only for elements over the main diagonal.


%% Input checks

% number of input args
if nargin ~= 2
    error('Function pruningFunction requires input arg "pValuesMatrix" and "matrixToPrune"!');
end

% check the size of input matrices
if size(pValuesMatrix, 1) ~= size(matrixToPrune, 1) || size(pValuesMatrix, 2) ~= size(matrixToPrune, 2)
    error('Function pruningFunction requires input marices "pValuesMatrix" and "matrixToPrune" of equal size!');
end


%% Order non-NaN matrix elements into a vector

% number of ROIs
roiNo = size(pValuesMatrix, 1);
% variable for storing ROI pairings we already calculated P value for
pastPairings = nan(roiNo*(roiNo-1)/2, 2);
% ROI pairing counter
counter = 0;
% vectorized p-values
vectorizedPValues = nan(roiNo, 1);

for roiOne = 1:roiNo  
    for roiTwo = 1:roiNo    
        % only consider P value if ROI numbers do not match and have not
        % been encountered before
        if roiOne ~= roiTwo && ~ismember([roiOne, roiTwo], pastPairings, 'rows') && ~ismember([roiTwo, roiOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [roiOne, roiTwo];
            
            vectorizedPValues(counter) = pValuesMatrix(roiOne, roiTwo);
            
        end  % if  
    end  % channelTwo
end  % channelOne


%% Set nonsignificant matrix elements into NaN

[h, pCrit] = fdr(vectorizedPValues);
prunedMatrix = matrixToPrune;
pastPairings = nan(roiNo*(roiNo-1)/2, 2);
counter = 0;

for roiOne = 1:roiNo
    for roiTwo = 1:roiNo
        % only consider P value if ROI numbers do not match and have not
        % been encountered before
        if roiOne ~= roiTwo && ~ismember([roiOne, roiTwo], pastPairings, 'rows') && ~ismember([roiTwo, roiOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [roiOne, roiTwo];
            if h(counter) ~= 1
                prunedMatrix(roiOne, roiTwo) = NaN;
            end
            
        end  % if    
    end  % channelTwo  
end  % channelOne

return
function [prunedMatrix] = pruningFunction(pValuesMatrix, matrixToPrune)

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
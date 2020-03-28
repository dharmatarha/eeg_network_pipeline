function pValues = surrConnStats(realConnData, surrConnData)

%% Calculate significance of connectivity values by comparing real to surrogate connectivity data
% 
% USAGE: pValues = surrConnStats(realConnData, surrConnData)
%
% Calculates empirical p-values (estimated significance) for each edge 
% (each entry in the connectivity matrix) by comparing real connectivity to
% connectivity in surrogate data.
%
% Inputs: 
% realConnData    - Numerical matrix (real), contains measured (real)
%                   connectivty data.
% surrConnData    - Numerical tensor (real), contains surrogate connectivity 
%                   connectivity data. Surrogate matrices are concatenated
%                   along the first dimension.
%
% Output:
% pValues         - Numerical matrix with the same size as input 
%                   matrix "realConnData", contains empirical p-values.
%
% Notes: It is assumed that the input matrix "realConnData" is symmetrical
% and p-values are calculated only for elements over the main diagonal.
%

%% Input checks

if ~ismatrix(realConnData) || ~isreal(realConnData)
    error('Function surrConnStats expects a numerical matrix of reals as input!');
end

if ~isreal(surrConnData)
    error('Function surrConnStats expects a numerical matrix of reals as input!');
end

% give a warning if there seems to be more
% variables than surrogates
if size(surrConnData, 2) > size(surrConnData, 1) || size(surrConnData, 3) > size(surrConnData, 1)
    warning(['There seems to be more variables than surrogates in ',...
        '- are you sure surrogate matrices are concatenated along the first',...
        'dimension of "surrConnData" ?'])
end


%% Calculation of p-values

% number of ROIs
roiNo = size(realConnData, 1);
% variable for storing ROI pairings we already calculated P value for
pastPairings = nan(roiNo*(roiNo-1)/2, 2);
% ROI pairing counter
counter = 0;
% variable for storing calculated p-values
pValues = nan(roiNo);

for roiOne = 1:roiNo
    
    for roiTwo = 1:roiNo
        
        % only calculate P value if ROI numbers do not match and have not
        % been encountered before
        if roiOne ~= roiTwo && ~ismember([roiOne, roiTwo], pastPairings, 'rows') && ~ismember([roiTwo, roiOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [roiOne, roiTwo];
            
            pValues(roiOne, roiTwo) = estimatedP(squeeze(realConnData(roiOne, roiTwo)), squeeze(surrConnData(:, roiOne, roiTwo)), 1);
            
        end  % if
        
    end  % channelTwo
    
end  % channelOne


return
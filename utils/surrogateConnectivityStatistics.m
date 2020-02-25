function pValues = surrogateConnectivityStatistics(measuredConnectivityData, surrogateConnectivityData)

%% Calculate surrogate connectivity statistics
% 
% USAGE: pValues = surrogateConnectivityStatistics(measuredConnectivityData, surrogateConnectivityData)
%
% Calculates empirical p-values of elements in the connectivity matrix 
% containing measured data based on the surrogate data.
%
% Works by calculating empirical p-values.
%
% Input: 
% measuredConnectivityData     - Numerical matrix (real), contains measured connectivty data.
% surrogateConnectivityData    - Numerical tensor (real), contains surrogate connectivity 
%                          data. Surrogate matrices are concatenated along the
%                          first dimension.
%
% Output:
% pValues                      - numerical matrix with the same size as input 
%                          matrix "measuredConnectivityData", contains empirical p-values.
%
% Notes: It is assumed that the input matrix "measuredConnectivityData" is symmetrical
% and p-values are calculated only for elements over the main diagonal.

%% Input checks

if ~ismatrix(measuredConnectivityData) || ~isreal(measuredConnectivityData)
    error('Function surrogateStatistics expects a numerical matrix of reals as input!');
end

if ~isreal(surrogateConnectivityData)
    error('Function surrogateStatistics expects a numerical matrix of reals as input!');
end

% give a warning if there seems to be more
% variables than surrogates
if size(surrogateConnectivityData, 2) > size(surrogateConnectivityData, 1) || size(surrogateConnectivityData, 3) > size(surrogateConnectivityData, 1)
    warning(['There seems to be more variables than surrogates in ',...
        '- are you sure surrogate matrices are concatenated along the first',...
        'dimension of "surrogateConnectivityData" ?'])
end


%% Calculation of p-values

% number of ROIs
roiNo = size(measuredConnectivityData, 1);
% variable for storing ROI pairings we already calculated P value for
pastPairings = nan(roiNo*(roiNo-1)/2, 2);
% ROI pairing counter
counter = 0;
% variable for storing calculated p-values
surrogatePs = nan(roiNo);

for roiOne = 1:roiNo
    
    for roiTwo = 1:roiNo
        
        % only calculate P value if ROI numbers do not match and have not
        % been encountered before
        if roiOne ~= roiTwo && ~ismember([roiOne, roiTwo], pastPairings, 'rows') && ~ismember([roiTwo, roiOne], pastPairings, 'rows')
            
            % remember pairing, adjust counter
            counter = counter+1;
            pastPairings(counter, :) = [roiOne, roiTwo];
            
            surrogatePs(roiOne, roiTwo) = estimatedP(squeeze(measuredConnectivityData(roiOne, roiTwo)), squeeze(surrogateConnectivityData(:, roiOne, roiTwo)), 1);
            
        end  % if
        
    end  % channelTwo
    
end  % channelOne

pValues = surrogatePs;

return
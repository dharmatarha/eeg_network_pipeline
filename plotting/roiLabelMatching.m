function [equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(roiLabelsData)
%% Function to check if a set of ROI labels match the standard ones expected for circleGraphPlotting
%
% USAGE: [equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(roiLabelsData)
%
% The function stores a long and a short version of the standard ROI label
% set expected and used by circleGraphPlotting. 
%
% If there is an exact match between the supplied and either of the 
% standard sets, the function returns equalFlag = true. If there is no 
% equality but there is a complete overlap, matchingSetsFlag is set to true. 
% In both cases, the function provides the matching standard label set in 
% roiLabelsPlotting. 
% 
% Input:
% roiLabelsData     - Cell array of char arrays, corresponding to the ROI
%           labels / names we want to match with expected ones.
%
% Outputs:
% equalFlag         - Boolean value, set to true if roiLabelsData is equal
%           to any of the two standard ROI label sets.
% matchingSetsFlag  - Boolean value, set to true if roiLabelsData overlaps
%           with any of the two standard ROI label sets.
% roiLabelsPlotting - If either equalFlag or matchingSetsFlag is true, it 
%           is a cell array of char arrays, containing the ROI label set 
%           corresponding to roilabelsData. Otherwise it is an empty array. 


%% Input checks

if nargin ~= 1
    error('Input arg "roiLabelsData" is required!');
end
if ~iscell(roiLabelsData) || ~isvector(roiLabelsData) || any(~cellfun(@ischar, roiLabelsData))
    error('Input arg "roiLabelsData" should be a cell vector of char arrays');
end
% transform to column vector if necessary
if ~iscolumn(roiLabelsData)
    roiLabelsData = roiLabelsData';
end


%% Standard ROI label sets for plotting (long and short versions)

% roi names used most often for our EEG datasets, grouped by lobules
roiLabelsPlottingLong = {' lateralorbitofrontal L ', ' medialorbitofrontal L ', ' parsorbitalis L ', ' parstriangularis L ', ' parsopercularis L ', ' rostralmiddlefrontal L ', ' caudalmiddlefrontal L ', ' superiorfrontal L ', ' precentral L ',...
    ' rostralanteriorcingulate L ', ' caudalanteriorcingulate L ', ' posteriorcingulate L ', ' isthmuscingulate L ',...
    ' transversetemporal L ', ' superiortemporal L ', ' middletemporal L ', ' inferiortemporal L ', ' entorhinal L ', ' parahippocampal L ', ' fusiform L ', ' insula L ',...
    ' supramarginal L ', ' inferiorparietal L ', ' superiorparietal L ', ' postcentral L ', ' paracentral L ', ' precuneus L ',...
    ' cuneus L ', ' lingual L ', ' pericalcarine L ', ' lateraloccipital L ',...
    ' lateraloccipital R ', ' pericalcarine R ', ' lingual R ', ' cuneus R ',...
    ' precuneus R ', ' paracentral R ', ' postcentral R ', ' superiorparietal R ', ' inferiorparietal R ', ' supramarginal R ',...
    ' insula R ', ' fusiform R ', ' parahippocampal R ', ' entorhinal R ', ' inferiortemporal R ', ' middletemporal R ', ' superiortemporal R ', ' transversetemporal R ',...
    ' isthmuscingulate R ', ' posteriorcingulate R ', ' caudalanteriorcingulate R ', ' rostralanteriorcingulate R ',...
    ' precentral R ', ' superiorfrontal R ', ' caudalmiddlefrontal R ', ' rostralmiddlefrontal R ', ' parsopercularis R ', ' parstriangularis R ', ' parsorbitalis R ', ' medialorbitofrontal R ', ' lateralorbitofrontal R ',...
    };
% shortened version
roiLabelsPlottingShort = {' latOrbFront L ', ' medOrbFront L ', ' parsOrb L ', ' parsTriang L ', ' parsOpercul L ', ' rostrMidFront L ', ' caudMidFront L ', ' supFront L ', ' precentral L ',...
    ' rostrAntCing L ', ' caudAntCing L ', ' postCing L ', ' isthmusCing L ',...
    ' transvTemp L ', ' supTemp L ', ' midTemp L ', ' infTemp L ', ' entorhinal L ', ' paraHippoc L ', ' fusiform L ', ' insula L ',...
    ' supraMarg L ', ' infPar L ', ' supPar L ', ' postcentral L ', ' paracentral L ', ' precuneus L ',...
    ' cuneus L ', ' lingual L ', ' periCalc L ', ' latOcc L ',...
    ' latOcc R ', ' periCalc R ', ' lingual R ', ' cuneus R ',...
    ' precuneus R ', ' paracentral R ', ' postcentral R ', ' supPar R ', ' infPar R ', ' supraMarg R ',...
    ' insula R ', ' fusiform R ', ' paraHippoc R ', ' entorhinal R ', ' infTemp R ', ' midTemp R ', ' supTemp R ', ' transvTemp R ',...
    ' isthmusCing R ', ' postCing R ', ' caudAntCing R ', ' rostrAntCing R ',...
    ' precentral R ', ' supFront R ', ' caudMidFront R ', ' rostrMidFront R ', ' parsOpercul R ', ' parsTriang R ', ' parsOrb R ', ' medOrbFront R ', ' latOrbFront R ',...
    };

% amount of shift needed for proper left-right arrangement in plot
shiftLabels = 15;
% shifting standard labels for plotting
roiLabelsPlottingLong = [roiLabelsPlottingLong(end-shiftLabels:end), roiLabelsPlottingLong(1:end-shiftLabels-1)]';
roiLabelsPlottingShort = [roiLabelsPlottingShort(end-shiftLabels:end), roiLabelsPlottingShort(1:end-shiftLabels-1)]';


%% Check if ROI label set from data match the sets used for plotting

% check for equality:
% set equalFlag if supplied arg "roiLabelsData" match either standard set of ROI labels
if isequal(roiLabelsPlottingLong , roiLabelsData)
    equalFlag = true;
    roiLabelsPlotting = roiLabelsPlottingLong;
elseif isequal(roiLabelsPlottingShort , roiLabelsData)
    equalFlag = true;
    roiLabelsPlotting = roiLabelsPlottingShort;
else
    equalFlag = false;
    roiLabelsPlotting = [];
end

% check for complete overlap:
% set matchingSetsFlag if supplied arg "roiLabelsData" overlaps completely with either standard set of ROI labels
if ~equalFlag  % only if sets were not equal
    intersectLabelsLong = intersect(roiLabelsData, roiLabelsPlottingLong);  % get intersect of supplied labels and standard ones
    intersectLabelsShort = intersect(roiLabelsData, roiLabelsPlottingShort);
    if isequal(intersectLabelsLong, sort(roiLabelsData))
        matchingSetsFlag = true; 
        roiLabelsPlotting = roiLabelsPlottingLong;
    elseif isequal(intersectLabelsShort, sort(roiLabelsData))
        matchingSetsFlag = true; 
        roiLabelsPlotting = roiLabelsPlottingShort;
    else
        matchingSetsFlag = false;
        roiLabelsPlotting = [];
    end
    
else
    matchingSetsFlag = true;
    
end  % if ~equalFlag


return
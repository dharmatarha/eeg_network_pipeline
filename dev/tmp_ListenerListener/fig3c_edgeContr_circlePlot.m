%% Script for Figure 3C
%
% The script generates the circle plot of the edges contributing the most
% to the within- vs. across-narratives effect.
%
% 

%% Base params

method = 'plv'; 
freq = 'alpha';

simMetric = 'corr';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load edge contributions data, determined by "method", "freq" and "thr"

edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_group_', thr, '.mat']);
data = load(edgeContrF);

% sanity check - were the edge contributions calculated with the requested
% similarity metric?
if ~strcmp(simMetric, data.simMetric)
    error(['Bad similarity metric! Edge contributions file specified simMetric: ', data.simMetric]);
end


%% Get cohen d values and p values in connectivity map format

ps = data.permRes.pEstAll; ps(ps==0) = 1 / data.permNo;
ds = data.permRes.cohendAll;

% in connectiviy map formats, by populating the upper triangles
pMap = nan(62); pMap(triu(true(62), 1)) = ps;
dMap = nan(62); dMap(triu(true(62), 1)) = ds;

% FDR
q = 0.05;
fdrMethod = 'bh';
[h, crit] = fdr(ps, q, fdrMethod);

% apply fdr crit to maps
pMapFdr = pMap; pMapFdr(pMap > crit) = 0;
dMapFdr = dMap; dMapFdr(pMap > crit) = 0;


%% Load ROI label names

labels = load('/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;  % short names


%% Load connectivity data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['group_surrResults_', freq, '_', method, '.mat']);
res = load(fileP);


%% Get mean connectivity, depending on "thr"

switch thr
    case 'unthr'
        connData = res.meanConn;
    case 'thrSub'
        % get the data from the averaged, subject-level thresholded matrices, both
        % negative and positive diffs
        connData = res.meanMaskedConnPos + res.meanMaskedConnNeg;   
    case 'thrGroup'
        connData = res.groupMaskedConnPos + res.groupMaskedConnNeg;
end  % switch
        
% get averaged conn matrix
meanConn = squeeze(mean(mean(connData, 1), 2));


%% Apply FDR

meanConnFdr = meanConn; 
meanConnFdr(pMap > crit) = 0;


%% Get new labels for reordering

newLabelsShort = {'latOrbFront L', 'medOrbFront L', 'parsOrb L', 'parsTriang L', 'parsOpercul L', 'rostrMidFront L', 'caudMidFront L', 'supFront L', 'precentral L',...
    'rostrAntCing L', 'caudAntCing L', 'postCing L', 'isthmusCing L',...
    'transvTemp L', 'supTemp L', 'midTemp L', 'infTemp L', 'entorhinal L', 'paraHippoc L', 'fusiform L', 'insula L',...
    'supraMarg L', 'infPar L', 'supPar L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'cuneus L', 'lingual L', 'periCalc L', 'latOcc L',...
    'latOcc R', 'periCalc R', 'lingual R', 'cuneus R',...
    'precuneus R', 'paracentral R', 'postcentral R', 'supPar R', 'infPar R', 'supraMarg R',...
    'insula R', 'fusiform R', 'paraHippoc R', 'entorhinal R', 'infTemp R', 'midTemp R', 'supTemp R', 'transvTemp R',...
    'isthmusCing R', 'postCing R', 'caudAntCing R', 'rostrAntCing R',...
    'precentral R', 'supFront R', 'caudMidFront R', 'rostrMidFront R', 'parsOpercul R', 'parsTriang R', 'parsOrb R', 'medOrbFront R', 'latOrbFront R',...
    };


%% Reorder connectivity and edge contributions matrices according to new labels

% amount of shift needed for proper left-right arrangement in plot
shiftLabels = 15;
newLabelsShort = [newLabelsShort(end-shiftLabels:end), newLabelsShort(1:end-shiftLabels-1)];

connMatrix = meanConnFdr;
[connMatrix, old2new] = matrixReorder(connMatrix, labels, newLabelsShort);
edgeWeights = dMapFdr(old2new, old2new);

% NaN values for zeros
connMatrix(connMatrix == 0) = nan;
edgeWeights(edgeWeights == 0) = nan;

% sanity check
if ~isequal(isnan(connMatrix), isnan(edgeWeights))
    error('Oops, non-matching NaN locations in connectivity and edge weight matrices!');
end


%% Circle plot

colorMap = 'jet';
mainFig = circleGraphPlot_edgeColorWeights(connMatrix, edgeWeights,... 
                                            colorMap, 0, newLabelsShort);


















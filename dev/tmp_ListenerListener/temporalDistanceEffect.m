%% Script for testing the effect of temporal distance between epochs
% The script loads the group-averaged connectivity similarity matrix 
% and tests if epoch distance has an effect on average similarity
%
%

%% Base params

method = 'ciplv'; 
freq = 'alpha';
simMethod = 'corr';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['mainCmpRes_', freq, '_', method, '.mat']);
res = load(fileP);

% sanity check, correct similarity method?
if ~strcmp(res.simMethod, simMethod)
    warning(['WRONG SIMILARITY METHOD? LOADED RESULTS FILE WAS OBTAINED WITH ', upper(res.simMethod), '!']);
end


%% Get connectivity similarity

connSim = res.(['groupRes_v2_', thr]).connSim;


%% Get distances for a given story

% correct distances will be in the lower triangle of the "dist" matrix
dist = nan(40);
for i = 1:39
    dist(1+i:41:end) = i;
end
% make matrix symmetrix
dist = tril(dist) + tril(dist)';

% vector format, only from upper triangle - linearization
distV = dist(triu(true(40), 1));


%% Extract similarity values per story

sims = nan(4, 40, 40);
for i = 1:4
    sims(i, :, :) = connSim((i-1)*40+1:i*40, (i-1)*40+1:i*40);
end
% vectorized version, from upper triangle
simsV = nan(4, 780);
for i = 1:4
    tmp = squeeze(sims(i, :, :));
    simsV(i, :) = tmp(triu(true(40), 1));
end


%% Get anovan model

% data in one vector
y = [simsV(1, :), simsV(2, :), simsV(3, :), simsV(4, :)]';

% grouping vars
allDist = repmat(distV, [4, 1]);
story = [repmat({'A'}, [780, 1]); repmat({'B'}, [780, 1]); repmat({'C'}, [780, 1]); repmat({'D'}, [780, 1])];
grouping = {allDist, story};

[P, T, STATS, TERMS] = anovan(y, grouping, 'continuous', [1], 'varnames', {'Epoch distance', 'Context'}, 'model', 'interaction');


%% Get GLM and LME

dataTable = table(allDist, story, y);

% Simple LM, with only the distance as predictor
LM = fitlm(dataTable, 'y~1+allDist', 'CategoricalVars', [false true false]);
% LM with stories as fixed predictors, with interaction
LM = fitlm(dataTable, 'y~1+allDist*story', 'CategoricalVars', [false true false]);
% LM with stories as fixed predictors, without interaction
LM = fitlm(dataTable, 'y~1+allDist+story', 'CategoricalVars', [false true false]);

% LME model, random intercepts and slopes for stories
LME = fitlme(dataTable, 'y~allDist+(allDist|story)', 'FitMethod', 'REML', 'CheckHessian', true, 'Verbose', true);
% LME model, random intercepts per stories
LME = fitlme(dataTable, 'y~allDist+(1|story)', 'FitMethod', 'REML', 'CheckHessian', true, 'Verbose', true);



%% Get LME/LM models without extreme dist values with small samples

adjTable = dataTable;
thresh = 25;
adjTable(adjTable.allDist>thresh, :) = [];

% Simple LM, with only the distance as predictor
LM = fitlm(adjTable, 'y~1+allDist', 'CategoricalVars', [false true false]);
% LM with stories as fixed predictors, with interaction
LM = fitlm(adjTable, 'y~1+allDist*story', 'CategoricalVars', [false true false]);
% LM with stories as fixed predictors, without interaction
LM = fitlm(adjTable, 'y~1+allDist+story', 'CategoricalVars', [false true false]);

% LME model, random intercepts and slopes for stories
LME = fitlme(adjTable, 'y~allDist+(allDist|story)', 'FitMethod', 'REML', 'CheckHessian', true, 'Verbose', true);
% LME model, random intercepts per stories
LME = fitlme(adjTable, 'y~allDist+(1|story)', 'FitMethod', 'REML', 'CheckHessian', true, 'Verbose', true);




%% Get results for only the second half

lastEpochs = 20;




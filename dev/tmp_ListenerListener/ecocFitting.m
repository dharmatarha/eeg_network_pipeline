function ecocFitting(freq)

%% ECOC script
% Script to calculate the "ideal" set of edges detecting a difference in
% context, using ECOC


%% Base params

method = 'ciplv'; 
% freq = 'alpha';
onlyPos = true;

simMetric = 'corr';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'groupv2thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load edge contributions data, determined by "method", "freq" and "thr"

if onlyPos
    edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '_onlyPos.mat']);
else
    edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '.mat']);
end
contrData = load(edgeContrF);

% sanity check - were the edge contributions calculated with the requested
% similarity metric?
if ~strcmp(simMetric, contrData.simMetric)
    error(['Bad similarity metric! Edge contributions file specified simMetric: ', contrData.simMetric]);
end


%% Get cohen d values and p values in connectivity map format

ps = contrData.permRes.pEstAll; ps(ps==0) = 1 / contrData.permNo;
ds = contrData.permRes.cohendAll;

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

% apply fdr crit to d / p values in vector format as well
dsFdr = ds;
dsFdr(ps > crit) = 0;
psFdr = ps;
psFdr(ps > crit) = 0;


%% Load ROI label names

labels = load('/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;  % short names


%% Load connectivity data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['group_surrResultsv2_', freq, '_', method, '.mat']);
res = load(fileP);


%% Get mean connectivity, depending on "thr"

switch thr
    case 'unthr'
        connData = res.meanConn;
    case 'groupv2thrSub'
        % get the data from the averaged, subject-level thresholded matrices, both
        % negative and positive diffs
        if onlyPos
            connData = res.meanMaskedConnPos;   
        else
            connData = res.meanMaskedConnPos + res.meanMaskedConnNeg; 
        end
    case 'groupv2thrGroup'
        if onlyPos
            connData = res.groupMaskedConnPos;
        else
            connData = res.groupMaskedConnPos + res.groupMaskedConnNeg;
        end
end  % switch
        
% get averaged conn matrix

meanConn = squeeze(mean(mean(connData, 1), 2));


%% Apply FDR to mean connectivity as well

meanConnFdr = meanConn; 
meanConnFdr(pMap > crit) = 0;


%% Fit ECOC

storyLabels = [repmat({'Story1'}, [40, 1]); repmat({'Story2'}, [40, 1]); repmat({'Story3'}, [40, 1]); repmat({'Story4'}, [40, 1])];
connData3d = cat(1, squeeze(connData(1,:,:, :)), squeeze(connData(2,:,:, :)), squeeze(connData(3,:,:, :)), squeeze(connData(4,:,:, :)));
%connDataLin = permute(connData3d, [2 3 1]); connDataLin = linearizeTrius(connDataLin, 1); connDataLin = connDataLin';
dsFdrSorted = sort(dsFdr, 'descend');
maxRuns = 100;
maxEdges = 300;

% var holding classification error rate results for each run and edge
% number
cer1 = zeros(maxRuns, maxEdges);

% ecoc params
kfold = 5;
codingType = 'onevsone';
learners = 'svm';
verbosity = 0;
opts = statset;
opts.UseParallel = false;

parfor topN = 1:maxEdges
    
        % select the top N edges as predictors
        topCritValue = dsFdrSorted(topN);
        ecocPred = connData3d(:, dMapFdr >= topCritValue);

        % temp var for collecting results in parfor
        tmpParfor = zeros(maxRuns, 1);
        
        for runN = 1:maxRuns
        
            % select the top N edges as predictors
            ecocObj = fitcecoc(ecocPred, storyLabels, 'Learners', learners,... 
                'Coding', codingType, 'KFold', kfold, 'Verbose', verbosity,...
                'Options', opts);
            
            tmpParfor(runN) = kfoldLoss(ecocObj);
        
        end  % for runNo
        
        % slice into results var
        cer1(:, topN) = tmpParfor;
        
end  % parfor topNo


saveF = ['ecocRes_', freq, '_', method];
save(saveF, 'cer1', 'freq', 'method');


end

% %% Different ECOC params
% 
% maxRuns = 100;
% maxEdges = 300;
% 
% % var holding classification error rate results for each run and edge
% % number
% cer2 = zeros(maxRuns, maxEdges);
% 
% % ecoc params
% kfold = 5;
% codingType = 'ternarycomplete';
% learners = 'svm';
% verbosity = 0;
% opts = statset;
% opts.UseParallel = false;
% 
% parfor topN = 1:maxEdges
%     
%         % select the top N edges as predictors
%         topCritValue = dsFdrSorted(topN);
%         ecocPred = connData3d(:, dMapFdr >= topCritValue);
% 
%         % temp var for collecting results in parfor
%         tmpParfor = zeros(maxRuns, 1);
%         
%         for runN = 1:maxRuns
%         
%             % select the top N edges as predictors
%             ecocObj = fitcecoc(ecocPred, storyLabels, 'Learners', learners,... 
%                 'Coding', codingType, 'KFold', kfold, 'Verbose', verbosity,...
%                 'Options', opts);
%             
%             tmpParfor(runN) = kfoldLoss(ecocObj);
%         
%         end  % for runNo
%         
%         % slice into results var
%         cer2(:, topN) = tmpParfor;
%         
% end  % parfor topNo
% 
% 
% 
%% Plot main results so far

m1 = mean(cer1); s1 = std(cer1);
% m2 = mean(cer2); s2 = std(cer2);
figure;
plot(m1, 'b-'); hold on; 
plot(m1+s1, 'b--'); plot(m1-s1, 'b--'); 
% plot(m2, 'r-'); plot(m2+s2, 'r--'); plot(m2-s2, 'r--'); 
hold off;



%% Filter the mean CER values and find the first local minima

filtOrder = 9;
cer1_medfilt = medfilt1(m1, filtOrder);
figure; plot(cer1_medfilt);
localMinima = find(diff(cer1_medfilt)>0);
disp([char(10), 'First few local minimas are with the first ', num2str(localMinima(1:5)), ' edges']);

finalN = localMinima(1); 
disp([char(10), 'Selected number of edges: ', num2str(finalN)]);

% prediction accuracies / errors can veiwed with kfoldPredict



%% Get connectivity matrix and edge list for topN edges selected in previous steps

finalCritValue = dsFdrSorted(finalN);
finalMap = dMapFdr >= finalCritValue;
finalDMap = dMapFdr; finalDMap(~finalMap) = 0;
finalConnMap = meanConn; finalConnMap(~finalMap) = 0;

% correlation between ds and connectivity strength values in final network
d_conn_corr = corr(finalDMap(finalMap), finalConnMap(finalMap));

% get matlab graph object
G = graph(finalConnMap, labels, 'upper');
% some useful metrics
table(G.centrality('degree'), G.centrality('betweenness'), labels')

% MST
G.minspantree('Root', 59).Edges


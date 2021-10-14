% Script for generating an empirical sample for each edge from the
% delete-one jackknifing estimates


%% Basic params

freq = 'alpha';
method = 'ciplv';
onlyPos = true;

simMetric = 'corr';

% type of thresholding (if any)
thr = 'groupv2thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';

if ~strcmp(freq, 'delta')
    subNo = 26;
else
    subNo = 25;
end

roiNo = 62;
edgeNo = roiNo*(roiNo-1)/2;


%% Load full-sample edge contributions data, determined by "method", "freq" and "thr"

if onlyPos
    edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '_onlyPos.mat']);
else
    edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '.mat']);
end
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
pMap = nan(roiNo); pMap(triu(true(roiNo), 1)) = ps;
dMap = nan(roiNo); dMap(triu(true(roiNo), 1)) = ds;

% FDR
q = 0.05;
fdrMethod = 'bh';
[h, crit] = fdr(ps, q, fdrMethod);

% apply fdr crit to maps
hMap = pMap > crit;
pMapFdr = pMap; pMapFdr(hMap) = 0;
dMapFdr = dMap; dMapFdr(hMap) = 0;


%% Collect jackknifing estimates

% preallocate
jps = zeros(edgeNo, subNo);
jds = jps;
jpMap = nan([roiNo, roiNo, subNo]);
jdMap = jpMap;


% go through delete-one jackknifed sample results
for subIdx = 1:subNo
    
    disp([char(10), 'Delete-one jackknife with subject ', num2str(subIdx), ' removed...']);
    
    % load from file
    if onlyPos
        edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '_onlyPos_jackknife', num2str(subIdx), '.mat']);
    else
        edgeContrF = fullfile(baseDir, 'edgeContr', ['edgeContr_', freq, '_', method, '_', thr, '.mat']);
    end
    jdata = load(edgeContrF);

    % sanity check - were the edge contributions calculated with the requested
    % similarity metric?
    if ~strcmp(simMetric, jdata.simMetric)
        error(['Bad similarity metric! Edge contributions file specified simMetric: ', data.simMetric]);
    end
    
    % collect p and d values
    tmp_ps = jdata.permRes.pEstAll; tmp_ps(tmp_ps==0) = 1 / jdata.permNo;  
    jps(:, subIdx) = tmp_ps;
    jds(:, subIdx) = jdata.permRes.cohendAll;    
    
    % in connectivity map format as well
    tmp_pmap = nan(roiNo); tmp_pmap(triu(true(roiNo), 1)) = jps(:, subIdx);
    tmp_dmap = nan(roiNo); tmp_dmap(triu(true(roiNo), 1)) = jds(:, subIdx);
    jpMap(:, :, subIdx) = tmp_pmap;
    jdMap(:, :, subIdx) = tmp_dmap;
    
end


%% Test if effect size samples are different from zero
% z-test for each edge
    
zh = nan(edgeNo, 1);
zp = zh;
zalpha = 0.05;

for i = 1:edgeNo
    
    [zh(i), zp(i)] = ztest(jds(i, :), 0, std(jds(i, :)), 'alpha', zalpha);
    
end
    
% Check if there were any edges for which the jackknifing (z-test) shows 
% no significant difference from zero, however, there was an overall effect
% on the level of the whole group
edgesToTrim = h & ~zh';

disp([char(10), 'No. of edges to trim as a result of jackknifing: ', num2str(sum(edgesToTrim))]);

% save results
saveF = fullfile(baseDir, ['jkknifeRes_edgeContr_', freq, '_', method, '_', thr, '_onlyPos.mat']);
save(saveF);
    
    
    
    
    
    
    
    


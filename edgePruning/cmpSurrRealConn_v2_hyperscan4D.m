function cmpSurrRealConn_v2_hyperscan4D(freq, varargin)
%% Function to compare actual (real) connectivity values to surrogate connectivity values
% Wrapper to be used on the hyperscan listener-listener (4D) data set.
%
% USAGE: cmpSurrRealConn_hyperscan4D(freq, 
%                                   dirName=pwd, 
%                                   subjects={'s02', 's03',...}, 
%                                   method='iplv')
%
% v2 : Version treating each epoch as an independent sample for the
% subject, thresholding based on the whole data.
%
% The function compares the actual (real) connectivity values of a subject
% to the surrogate values. 
% Real connectivity values are - most likely - from the outcome of 
% "connectivityWrapperReal_hyperscan4D.m", that is, from subject-level 
% .mat files conforming to the file name pattern
% SUBJECTID_FREQUENCYBAND_METHOD.mat (e.g. "s03_alpha_iplv.mat").
% Surrogate values are expected to be the outcomes from
% "surrEdgeEstimationReal_hyperscan4D.m", that is, stored in files
% conforming to the naming pattern
% SUBJECTID_FREQUENCYBAND_surrEdgeEstReal_METHOD.mat. 
%
% NOTE that real connectivity data should be stored in 4D array "connRes",
% while surrogate data file should contain 5D arrays "surrNormalMu",
% "surrNormalH", "surrNormalP", "surrNormalSigma", and cell arrays 
% "truncated" and "method".
%
% Logic:
% (1) Loop through subjects
% (2) For each subject, load actual and surrogate connectivity data
% (3) Check if surrogate data could be fitted well enough with (truncated)
% normals during surrogate calculations
% (4) For each edge, compare the mean of the actual values to the mean of 
% the surrogate values, and derive a significance value
% (5) Perform FDR on the resulting connectivity map.
%
% The outputs are saved out into a file
% 'group_surrResultsv2_FREQUENCYBAND_METHOD.mat'.
%
% Mandatory inputs:
% freq              - Char array, one of {'delta', 'theta', 'alpha', 
%                   'beta', 'gamma'}. Frequency band which is reflected
%                   in surrogate connectivity file names.
%
% Optional inputs:
% dirName           - Char array, path to folder containing the folder of 
%                   frequency band data. E.g. subject-level data files if 
%                   "freq" == 'alpha' are expected to be located under 
%                   [dirName, '/', freq, '/']. Default is current working
%                   directory (pwd).
% subjects          - Cell array, list of subjects whose data we process. 
%                   Each cell contains a subject ID also used in the 
%                   filenames (e.g. 's01' for files like 
%                   's01_alpha_plv.mat'). The default is the subject list
%                   for the listener-listener hyperscan data set:
%                   subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%                   ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%                   's21','s22','s23','s24','s25','s26','s27','s28'}
% method            - Char array, one of {'plv', 'iplv', 'ciplv', 
%                   'ampCorr', 'orthAmpCorr'}. Connectivity method,
%                   reflected in connectivity file names. Defaults to
%                   'iplv'. 
%
% Outputs:
%
% The output is saved into a file named 
% "group_surrResultsv2_FREQUENCYBAND_METHOD.mat".
% It contains the following variables:
%
% realConn
%
% realConnP             - 4D numeric array, contains the probabilities 
%                       of the real connectivity values based on the normal
%                       distributions of the surrogate values. Sized 
%                       [subjects X epochs X rois X rois].
%
% maskedConnPos         - 4D numeric array, connectivity values after
%                       thresholding based on surrogate connectivity
%                       values.
%
% maskedConnNeg         - 4D numeric array, connectivity values after
%                       thresholding based on surrogate connectivity
%                       values.
%
% meanConn
%
% meanMaskedConnPos
%
% meanMaskedConnNeg
%
% diffDirection
%
% criticalP             - 2D numeric array, contains the critical p values
%                       for each epoch from FDR (q = .05), used for 
%                       thresholding. Sized [subjects X epochs].
%
% survivalRate          - 2D numeric array, contains the ratio of edges
%                       surviving the FDR-based thresholding in each epoch.
%                       Sized [subjects X epochs].
%
% failedFitRate
%
%
% surrNormalMu          - 4D numeric array, contains the "mu"s of the
%                       normal distributions fitted to the surrogate 
%                       connectivity values. Sized [subjects X epochs X
%                       rois X rois].
% surrNormalSigma       - 4D numeric array, contains the "sigma"s of the
%                       normal distributions fitted to the surrogate 
%                       connectivity values. Sized [subjects X epochs X
%                       rois X rois].
% surrNormalP           - 4D numeric array, contains the probabilities 
%                       from the Kolmogorov-Smirnov tests on the normal 
%                       fits of the surrogate connectivity values. Sized 
%                       [subjects X epochs X rois X rois].
%
% groupP
%
% groupCritP 
%
% groupDiffDir 
% 
% groupMaskedConnPos
%
% groupMaskedConnNeg
%
% groupSurvRate
%


%% Input checks

% check no. of inputs
if ~ismember(nargin, 1:4)
    error('Function cmpSurrRealConn_hyperscan4D requires input arg "freq", while "dirName", "subjects" and "method" are optional!');
end
% check mandatory input
if ~ismember(freq, {'alpha', 'beta', 'gamma', 'delta', 'theta'})
    error('Input arg "freq" should be one of {"alpha", "beta", "gamma", "delta", "theta"}!');
end
% check optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ~exist('dirName', 'var') && exist(varargin{v}, 'dir')
            dirName = varargin{v};
        elseif iscell(varargin{v}) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'ciplv', 'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}) && ~exist('method', 'var')
            method = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "dirName", "subjects" or "method"!');
        end
    end
end
% assign defaults
if ~exist('dirName', 'var')
    dirName = pwd;
end
if ~exist('subjects', 'var')
  subjects = {'s02','s03','s04','s05','s06','s07','s08','s09',...
      's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
      's21','s22','s23','s24','s25','s26','s27','s28'};
end
if ~exist('method', 'var')
    method = 'iplv';
end 

% user message
disp([char(10), 'Called function cmpSurrRealConn_v2_hyperscan4D with inputs: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Data folder: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Subject list: ']);
disp(subjects);


%% Settings, hardcoded params

subNo = length(subjects);
truncateBounds = [0 1];  % bounds if the surrogate prob distributions were truncated
fdrQ = 0.05;  % q for FDR
fdrType = 'bh';  % type of FDR correction to apply
surrDimPerm = [4, 3, 1, 2];  % vector for permuting the dimensions of the surrogate data so it matches the dimensions of the real connectivity data, from ROIs X ROIs X epochs X conditions to conditions X epochs X ROIs X ROIs
saveFile = fullfile(dirName, freq, ['group_surrResultsv2_', freq, '_', method, '.mat']);  % path for saving out results


%% Check size of data in first file

% char subject ID
subID = subjects{1};
% expected real connectivity data file name
subRealFile = fullfile(dirName, freq, [subID, '_', freq, '_', method, '.mat']);
if ~exist(subRealFile, 'file')
    error(['Cannot find connectivity data for subject ', subID, ' at ', subRealFile, '!']);
end
% load real connectivity data
connData = load(subRealFile);
% store the size of connectivity data 
if numel(size(connData.connRes)) ~= 4
    error(['Expected 4D connectivity array in file ', subRealFile, '!']);
end
[condNo, epochNo, roiNo, ~] = size(connData.connRes);


%% Preallocate output vars

% vector holding the ratios of failed (truncated) normal fits on surrogate
% data
failedFitRate = nan(subNo, 1);
% array of significance values for actual connectivity values
realConnP = nan([subNo, roiNo, roiNo]);
% array for the direction fo difference between real and surrogate
% connectivity values
diffDirection = realConnP;
% array for the FDR corrected, significantly stronger edges
maskedConnPos = nan([subNo, condNo, epochNo, roiNo, roiNo]);
% array for the FDR corrected, significantly weaker edges
maskedConnNeg = nan([subNo, condNo, epochNo, roiNo, roiNo]);
% array for collecting all subjects' raw connectivity values
realConn = nan([subNo, condNo, epochNo, roiNo, roiNo]);
% array of critical P values
criticalP = nan([subNo, 1]);
% array of FDR-surviving edge rates
survivalRate = criticalP;
% array for subject-averaged normal distribution mu, sigma and connectivity
surrMuSubAvg = nan([subNo, roiNo, roiNo]);
surrSigmaSubAvg = nan([subNo, roiNo, roiNo]);
connDataSubAvg = nan([subNo, roiNo, roiNo]);

% user message
disp([char(10), 'Prepared everything, starting loop across subjects.']);
    

%% Loop across subjects

for subIdx = 1:subNo
    
    % subject-level clock
    subClock = tic;
    
    % char subject ID
    subID = subjects{subIdx};
    
    % user message
    disp([char(10), 'Working on data from subject ', subID, '...']);
    
    
    %% Load real connectivity data
    
    % expected real connectivity data file name
    subRealFile = fullfile(dirName, freq, [subID, '_', freq, '_', method, '.mat']);
    if ~exist(subRealFile, 'file')
        error(['Cannot find connectivity data for subject ', subID, ' at ', subRealFile, '!']);
    end
    
    % load real connectivity data
    connData = load(subRealFile);
    % sanity checks
    if ~strcmp(connData.method, method)
        error('Connectivity measure in real connectivity data file does not match the measure specified in input args!');
    end
    if ~strcmp(connData.freq, freq)
        error('Frequency band in real connectivity data file does not match the band specified in input args!');
    end
    
    % store the size of connectivity data for the first subject, test for
    % others
    if subIdx == 1
        if numel(size(connData.connRes)) ~= 4
            error(['Expected 4D connectivity array in file ', subRealFile, '!']);
        else
            [condNo, epochNo, roiNo, ~] = size(connData.connRes);
        end
    else 
        if ~isequal([condNo, epochNo, roiNo, roiNo], size(connData.connRes))
            error('Connectivity data has unexpected size, investigate!');
        end
    end
    
    % keep only the connectivity data part
    connData = connData.connRes;  % dims conditions X epochs X ROIs X ROIs
    
    % user message
    disp('Loaded actual connectivity data');
    
    
    %% Load surrogate connectivity data
        
    % expected surrogate connectivity data file name
    subSurrFile = dir(fullfile(dirName, freq, [subID, '_', freq, '_surrEdgeEstReal*_', method, '*.mat']));
    if numel(subSurrFile) ~= 1
        error(['There are none or too many surrogate data file(s) for subject ', subID, '!']);
    end
    
    % load surrogate connectivity data
    surrData = load(fullfile(subSurrFile.folder, subSurrFile.name));   
    % sanity check
    if ~strcmp(surrData.subject, subID)
        error('Subject ID in surrogate connectivity data file does not match the ID of current subject!');
    end

    % if there are surrogate values for multiple methods, extract the
    % index corresponding to the current method
    if numel(surrData.method) ~= 1 && ismember(method, surrData.method)
        methodIdx = find(strcmp(method, surrData.method));
    elseif numel(surrData.method) == 1 && ismember(method, surrData.method)  % if there is only one method and it matches input arg "method"
        methodIdx = 1;
    else
        error(['Cannot match input arg "method" to the methods in surrogate data file at ', subSurrFile, '!']);
    end
        
    % was the surrogate data fitted with truncated normals? 
    truncatedFlag = surrData.truncated(methodIdx);
    if ~truncatedFlag && ismember(method, {'plv', 'iplv' ,'ciplv'})
        warning('Method is phase based but surrogate normals are not truncated (truncatedFlag is false)!');
    end
    
    % extract surrogate values we will use  
    surrH = squeeze(surrData.surrNormalH(methodIdx, :, :, :, :));  % dims ROIs X ROIs X epochs X conditions
    surrMu = squeeze(surrData.surrNormalMu(methodIdx, :, :, :, :));  % dims ROIs X ROIs X epochs X conditions
    surrSigma = squeeze(surrData.surrNormalSigma(methodIdx, :, :, :, :));  % dims ROIs X ROIs X epochs X conditions

    % get rate of successful fits of surrogate data
    edgesNo = condNo*epochNo*roiNo*(roiNo-1)/2;  % only upper triangles should hold values
    failedFitRate(subIdx) = sum(surrH(:), 'omitnan')/edgesNo;
    
    % user message
    disp('Loaded surrogate connectivity data');
    disp(['Ratio of failed fits on surrogate data: ',... 
        num2str(100*failedFitRate(subIdx)), '%']);
    

    %% Prepare connectivity and surrogate data for groupSurrStats
    % Reshape them into 3D arrays with dims ROIs X ROIs X epochs
    
    connDataReshaped = permute(connData, [3 4 2 1]);  % dims ROIs X ROIs X epochs X conditions
    connDataReshaped = reshape(connDataReshaped, [roiNo, roiNo, epochNo*condNo]);  % dims ROIs X ROIs X epochs*conditinos
    surrMuReshaped = reshape(surrMu, [roiNo, roiNo, epochNo*condNo]);  % dims ROIs X ROIs X epochs*conditinos
    surrSigmaReshaped = reshape(surrSigma, [roiNo, roiNo, epochNo*condNo]);  % dims ROIs X ROIs X epochs*conditinos


    %% Get signifance values as derived from the combination of all epochs
    
    if truncatedFlag
        [p, d] = groupSurrStats(connDataReshaped, surrMuReshaped, surrSigmaReshaped, truncateBounds, 'silent');
    else
        [p, d] = groupSurrStats(connDataReshaped, surrMuReshaped, surrSigmaReshaped, 'silent');
    end
    
    % store in aggregate arrays
    realConnP(subIdx, :, :) = p;
    diffDirection(subIdx, :, :) = d;
    
    % FDR 
    realPs = p(:); realPs(isnan(realPs)) = [];
    [~, criticalP(subIdx)] = fdr(realPs, fdrQ, fdrType);

    % create masked connectivity tensor
    pTensor = repmat(p, [1, 1, epochNo, condNo]);  % one layer of same p values for each epoch in each cond
    pTensor = permute(pTensor, [4 3 1 2]);  % dims conditions X epochs X ROIs X ROIs
    connDataMasked = connData;
    connDataMasked(pTensor > criticalP(subIdx)) = 0;
    
    % get positively different part
    dTensor = repmat(d, [1, 1, epochNo, condNo]);  % one layer of same p values for each epoch in each cond
    dTensor = permute(dTensor, [4 3 1 2]);  % dims conditions X epochs X ROIs X ROIs    
    connPos = connDataMasked;
    connPos(dTensor ~= 1) = 0;  % only edges with positive difference (real > surrogate)
    maskedConnPos(subIdx, :, :, :, :) = connPos;
    % get negatively different part
    connNeg = connDataMasked;
    connNeg(dTensor ~= -1) = 0;  % only edges with positive difference (real > surrogate)
    maskedConnNeg(subIdx, :, :, :, :) = connNeg;
    
    % get ratio of edges below the threshold (edges surviving pruning)
    edgesBelowCrit = p <= criticalP(subIdx);
    survivalRate(subIdx) = sum(edgesBelowCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);    
    
    % user message
    disp(['Ratio of edges surviving thresholding: ',... 
        num2str(survivalRate(subIdx)*100), '%']);
    
    
    %% Get average surrogate distributions for each subject (across all epochs)
    
    surrMuSubAvg(subIdx, :, :) = mean(surrMuReshaped, 3, 'omitnan');
    surrSigmaSubAvg(subIdx, :, :) = sqrt(sum(surrSigmaReshaped.^2, 3)./((epochNo*condNo)^2));
    connDataSubAvg(subIdx, :, :) = squeeze(mean(mean(connData, 1, 'omitnan'), 2, 'omitnan'));
    
    
    %% Collect group-level values
    
    realConn(subIdx, :, :, :, :) = connData;
    
    % report elapsed time
    subTime = round(toc(subClock), 3);
    disp(['Done with subject ', subID, '. Elapsed time: ', num2str(subTime), ' secs']);
    
    
end  % for subIdx

% get group averages
meanConn = squeeze(mean(realConn, 1));
meanMaskedConnPos = squeeze(mean(maskedConnPos, 1));
meanMaskedConnNeg = squeeze(mean(maskedConnNeg, 1));


%% Get also group-level thresholding results

% user message
disp([char(10), 'Finished with all subjects, calculating group-level thresholds now']);

% clock for group-level calculations
groupClock = tic;

% dimension rearrangements for calling groupSurrStats
tmp_conn = permute(connDataSubAvg, [2 3 1]);
tmp_mu = permute(surrMuSubAvg, [2 3 1]);
tmp_sigma = permute(surrSigmaSubAvg, [2 3 1]);

% Get p-values
if truncatedFlag
    [groupP, groupDiffDir] = groupSurrStats(tmp_conn, tmp_mu, tmp_sigma, truncateBounds, 'silent');
else
    [groupP, groupDiffDir] = groupSurrStats(tmp_conn, tmp_mu, tmp_sigma, 'silent');
end

% FDR
groupPs = groupP(:); groupPs(isnan(groupPs)) = [];
[~, groupCritP] = fdr(groupPs, fdrQ, fdrType);

% create masked connectivity tensor
pTensor = repmat(groupP, [1, 1, epochNo, condNo]);  % one layer of same p values for each epoch in each cond
pTensor = permute(pTensor, [4 3 1 2]);  % dims conditions X epochs X ROIs X ROIs
meanConnMasked = meanConn;
meanConnMasked(pTensor > groupCritP) = 0;     

% get positively different part
dTensor = repmat(groupDiffDir, [1, 1, epochNo, condNo]);  % one layer of same p values for each epoch in each cond
dTensor = permute(dTensor, [4 3 1 2]);  % dims conditions X epochs X ROIs X ROIs    
groupMaskedConnPos = meanConnMasked;
groupMaskedConnPos(dTensor ~= 1) = 0;  % only edges with positive difference (real > surrogate)

% get negatively different part
groupMaskedConnNeg = meanConnMasked;
groupMaskedConnNeg(dTensor ~= -1) = 0;  % only edges with positive difference (real > surrogate)

% get ratio of edges below the threshold (edges surviving pruning)
edgesBelowCrit = groupP <= groupCritP;
groupSurvRate = sum(edgesBelowCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);   
            
% get elapsed time
groupTime = round(toc(groupClock), 3);

% user message
disp(['Done. Group-level thresholding took ', num2str(groupTime), ' secs']);

            
%% Saving, cleaning up

save(saveFile, 'realConn', 'realConnP', 'maskedConnPos', 'maskedConnNeg',...
    'meanConn', 'meanMaskedConnPos', 'meanMaskedConnNeg',...
    'surrMuSubAvg', 'surrSigmaSubAvg', 'connDataSubAvg', ...
    'diffDirection', 'criticalP', 'survivalRate', 'failedFitRate',... 
    'groupP', 'groupCritP', 'groupDiffDir', 'groupMaskedConnPos',...
    'groupMaskedConnNeg', 'groupSurvRate');

% user message
disp([char(10), 'Saved out results to ', saveFile]);
disp('Done, finished, over, curtains.');


return
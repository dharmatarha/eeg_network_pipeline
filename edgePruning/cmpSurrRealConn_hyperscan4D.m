function cmpSurrRealConn_hyperscan4D(freq, varargin)
%% Function to compare actual (real) connectivity values to surrogate connectivity values
% Wrapper to be used on the hyperscan listener-listener (4D) data set.
%
% USAGE: cmpSurrRealConn_hyperscan4D(freq, 
%                                   dirName=pwd, 
%                                   subjects={'s02', 's03',...}, 
%                                   method='iplv')
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
% (4) Compare actual values to surrogate values and derive
% significance values
% (5) And perform FDR on the epoch-level
% (6) Do the same thresholding on group-level as well, for all epochs
%
% The outputs are saved out into a file
% 'group_surrResults_FREQUENCYBAND_METHOD.mat'.
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
% "group_surrResults_FREQUENCYBAND_METHOD.mat".
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
disp([char(10), 'Called function cmpSurrRealConn_hyperscan4D with inputs: ',...
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
surrDimPerm = [4, 3, 1, 2];  % vector for permuting the dimensions of the surrogate data so it matches the dimensions of the real connectivity data
saveFile = fullfile(dirName, freq, ['group_surrResults_', freq, '_', method, '.mat']);  % path for saving out results


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
realConnP = nan([subNo, condNo, epochNo, roiNo, roiNo]);
% array for the direction fo difference between real and surrogate
% connectivity values
diffDirection = realConnP;
% array for the FDR corrected, significantly stronger edges
maskedConnPos = realConnP;
% array for the FDR corrected, significantly weaker edges
maskedConnNeg = realConnP;
% array for collecting all subjects' raw connectivity values
realConn = realConnP;
% arrays for collecting all subjects' surrogate params
surrNormalH = realConnP;
surrNormalMu = realConnP;
surrNormalSigma = realConnP;
surrNormalP = realConnP;
% array of critical P values
criticalP = nan([subNo, condNo, epochNo]);
% array of FDR-surviving edge rates
survivalRate = criticalP;

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
    connData = connData.connRes;
    
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
    elseif surrData.method == 1 && ismember(method, surrData.method)  % if there is only one method and it matches inpuat arg "method"
        methodIdx = 1;
    else
        error(['Cannot match input arg "method" to the methods in surrogate data file at ', subSurrFile, '!']);
    end
        
    % was the surrogate data fitted with truncated normals? 
    truncatedFlag = surrData.truncated(methodIdx);
    if ~truncatedFlag && ismember(method, {'plv', 'iplv'})
        warning('Method is phase based but surrogate normals are not truncated (truncatedFlag is false)!');
    end
    
    % extract surrogate values we will use   
    surrH = squeeze(surrData.surrNormalH(methodIdx, :, :, :, :));
    surrMu = squeeze(surrData.surrNormalMu(methodIdx, :, :, :, :));
    surrSigma = squeeze(surrData.surrNormalSigma(methodIdx, :, :, :, :));
    surrP = squeeze(surrData.surrNormalP(methodIdx, :, :, :, :));

    % permute the order of dimensions on surrogate values
    surrH = permute(surrH, surrDimPerm);
    surrMu = permute(surrMu, surrDimPerm);
    surrSigma = permute(surrSigma, surrDimPerm);
    surrP = permute(surrP, surrDimPerm);
    
    % sanity check for size
    if ~isequal([condNo, epochNo, roiNo, roiNo], size(surrH))
        error('Surrogate values have unexpected size after dimension reordering, investigate!');
    end

    % get rate of successful fits of surrogate data
    edgesNo = condNo*epochNo*roiNo*(roiNo-1)/2;  % only upper triangles should hold values
    failedFitRate(subIdx) = sum(surrH(:), 'omitnan')/edgesNo;
    
    % user message
    disp('Loaded surrogate connectivity data');
    disp(['Ratio of failed fits on surrogate data: ',... 
        num2str(100*failedFitRate(subIdx)), '%']);
    
    
    %% Determine the direction of difference between real and surrogate connectivity data
    
    tmpDiff = nan([condNo, epochNo, roiNo, roiNo]);
    tmpDiff(connData<surrMu) = -1;
    tmpDiff(connData==surrMu) = 0;
    tmpDiff(connData>surrMu) = 1;
    diffDirection(subIdx, :, :, :, :) = tmpDiff;
    

    %% Get significance values of real connectivity values
    
    for condIdx = 1:condNo
        
        for epochIdx = 1:epochNo
            
            for roi1 = 1:roiNo
                for roi2 = 1:roiNo
                    % only for upper triangle of connetivity matrix
                    if roi1 < roi2

                        % if the surrogate normal was truncated, create a
                        % truncated probablity distribution for testing
                        % against
                        if truncatedFlag
                            
                            % get simple normal first with the surrogate
                            % params
                            pd = makedist('normal', 'mu', surrMu(condIdx, epochIdx, roi1, roi2),... 
                                'sigma', surrSigma(condIdx, epochIdx, roi1, roi2));  % returns a probability distribution object
                            
                            % truncate
                            pd = truncate(pd, truncateBounds(1), truncateBounds(2));
                            
                            % get p-value from the cumulative version
                            normalP = pd.cdf(connData(condIdx, epochIdx, roi1, roi2));

                        % Case of using simple normal distribution
                        elseif ~truncatedFlag

                            % get p-value based on normal distribution of surrogate
                            % values
                            normalP = normcdf(connData(condIdx, epochIdx, roi1, roi2),... 
                                                surrMu(condIdx, epochIdx, roi1, roi2),... 
                                                surrSigma(condIdx, epochIdx, roi1, roi2));

                        end  % if truncatedFlag

                        % convert for a 2-two-sided test, apply correction as
                        % well
                        if normalP > 0.5
                            normalP = 1-normalP;
                        end                                         
                        realConnP(subIdx, condIdx, epochIdx, roi1, roi2) = normalP*2;
                        
                    end  % if roi1 < roi2
                    
                end  % for roi2
            end  % for roi1 
            
            % FDR correction on epoch-level
            realPs = realConnP(subIdx, condIdx, epochIdx, :, :); 
            realPs = realPs(:); realPs(isnan(realPs)) = [];
            [~, criticalP(subIdx, condIdx, epochIdx)] = fdr(realPs, fdrQ, fdrType);             

            % create masked connectivity tensor
            epochConnData = squeeze(connData(condIdx, epochIdx, :, :));
            epochP = squeeze(realConnP(subIdx, condIdx, epochIdx, :, :));
            epochConnData(epochP > criticalP(subIdx, condIdx, epochIdx)) = 0;  % only edges surviving FDR
            epochPos = epochConnData; 
            epochPos(squeeze(diffDirection(subIdx, condIdx, epochIdx, :, :)) ~= 1) = 0;  % only edges with positive difference (real > surrogate)
            maskedConnPos(subIdx, condIdx, epochIdx, :, :) = epochPos;
            epochNeg = epochConnData; 
            epochNeg(squeeze(diffDirection(subIdx, condIdx, epochIdx, :, :)) ~= -1) = 0;  % only edges with positive difference (real > surrogate)
            maskedConnNeg(subIdx, condIdx, epochIdx, :, :) = epochNeg;            
            % get ratio of edges below the threshold (edges surviving the
            % pruning)
            edgesBelowCrit = epochP <= criticalP(subIdx, condIdx, epochIdx);
            survivalRate(subIdx, condIdx, epochIdx) = sum(edgesBelowCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);
             
        end  % for epochIdx
        
    end  % for condIdx


    %% Collect group-level values
    
    realConn(subIdx, :, :, :, :) = connData;
    surrNormalH(subIdx, :, :, :, :) = surrH;
    surrNormalMu(subIdx, :, :, :, :) = surrMu;
    surrNormalSigma(subIdx, :, :, :, :) = surrSigma;
    surrNormalP(subIdx, :, :, :, :) = surrP;
    
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
disp([char(10), 'Finished with all subjects, calculating group-level thresholds now for each epoch']);

% clock for group-level calculations
groupClock = tic;

groupP = nan([condNo, epochNo, roiNo, roiNo]);
groupMaskedConnPos = groupP;
groupMaskedConnNeg = groupP;
groupDiffDir = groupP;
groupCritP = nan([condNo, epochNo]);
groupSurvRate = groupCritP;

for condIdx = 1:condNo
    for epochIdx = 1:epochNo
        
        % epoch-level group data
        tmpConn = squeeze(realConn(:, condIdx, epochIdx, :, :));
        tmpSurrMu = squeeze(surrNormalMu(:, condIdx, epochIdx, :, :));
        tmpSurrSigma = squeeze(surrNormalSigma(:, condIdx, epochIdx, :, :));
        % permute dimensions for calling groupSurrStats
        tmpConn = permute(tmpConn, [2 3 1]);
        tmpSurrMu = permute(tmpSurrMu, [2 3 1]);
        tmpSurrSigma = permute(tmpSurrSigma, [2 3 1]);
        % get significance values for each edge in given epoch
        if ~truncatedFlag
            [groupP(condIdx, epochIdx, :, :), groupDiffDir(condIdx, epochIdx, :, :)] = groupSurrStats(tmpConn, tmpSurrMu, tmpSurrSigma, 'silent');
        else
            [groupP(condIdx, epochIdx, :, :), groupDiffDir(condIdx, epochIdx, :, :)] = groupSurrStats(tmpConn, tmpSurrMu, tmpSurrSigma, truncateBounds, 'silent');
        end
        
        % fdr correction per epoch
        tmpPs = groupP(condIdx, epochIdx, :, :); 
        tmpPs = tmpPs(:); tmpPs(isnan(tmpPs)) = [];
        [~, groupCritP(condIdx, epochIdx)] = fdr(tmpPs, fdrQ, fdrType); 
        
        % create masked connectivity tensor
        epochData = squeeze(meanConn(condIdx, epochIdx, :, :));
        epochP = squeeze(groupP(condIdx, epochIdx, :, :));
        epochData(epochP > groupCritP(condIdx, epochIdx)) = 0;  % only edges surviving FDR
        epochPos = epochData; 
        epochPos(squeeze(groupDiffDir(condIdx, epochIdx, :, :)) ~= 1) = 0;  % only edges with positive difference (real > surrogate)
        groupMaskedConnPos(condIdx, epochIdx, :, :) = epochPos;
        epochNeg = epochData; 
        epochNeg(squeeze(diffDirection(subIdx, condIdx, epochIdx, :, :)) ~= -1) = 0;  % only edges with positive difference (real > surrogate)
        groupMaskedConnNeg(condIdx, epochIdx, :, :) = epochNeg;            
        % get ratio of edges below the threshold (edges surviving the
        % pruning)
        edgesBelowCrit = epochP <= groupCritP(condIdx, epochIdx);
        groupSurvRate(condIdx, epochIdx) = sum(edgesBelowCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);        
        
    end  % for epochIdx
end  % for condIdx
            
% get elapsed time
groupTime = round(toc(groupClock), 3);

% user message
disp(['Done. Group-level thresholding took ', num2str(groupTime), ' secs']);

            
%% Saving, cleaning up

save(saveFile, 'realConn', 'realConnP', 'maskedConnPos', 'maskedConnNeg',...
    'meanConn', 'meanMaskedConnPos', 'meanMaskedConnNeg',...
    'diffDirection', 'criticalP', 'survivalRate', 'failedFitRate',... 
    'surrNormalH', 'surrNormalMu', 'surrNormalSigma', 'surrNormalP',...
    'groupP', 'groupCritP', 'groupDiffDir', 'groupMaskedConnPos',...
    'groupMaskedConnNeg', 'groupSurvRate');

% user message
disp([char(10), 'Saved out results to ', saveFile]);
disp('Done, finished, over, curtains.');


return


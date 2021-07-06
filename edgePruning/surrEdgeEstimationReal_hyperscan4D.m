function surrEdgeEstimationReal_hyperscan4D(freq, varargin)

%% Estimate distribution of edge weights (adj. matrix) for surrogate / null data
% Wrapper to be used on the hyperscan listener-listener (4D) data set.
%
% USAGE: surrEdgeEstimationReal_hyperscan4D(freq,
%                                           dirName = pwd, 
%                                           subjects = {'s02', 's03', ...}, 
%                                           method = 'iplv',
%                                           dRate = 1,
%                                           surrNo = 1000,
%                                           truncated = 'truncated',
%                                           failedFitAction='saveResults')
%
% Fits edge weights for phase-scrambling-based surrogate data with a 
% (truncated) normal distribution, estimating the parameters (mean and std). 
% USES PARFOR!
%
% Connectivity is measured by calling functions (plv, iplv, ciplv, ampCorr or orthAmpCorr)
% defined outside this script. 
%
% NOTE that for envelope correlation measures you should fit a
% non-truncated normal to the (Fisher-z tranformed) surrogate values.
%
% NOTE that for envelope correlation measures we perform a hard-coded
% lowpass filtering with 10 Hz passbandedge on the envelopes before calculating 
% correlations. This step - currently - requires Signal Processing
% Toolbox!!!
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'SUBJECTNUMBER_FREQUENCYBAND_surrEdgeEstReal_METHOD.mat'.
%
% Output files contain the parameters of the fitted (truncated) normal
% distributions and the results of goodness-of-fit tests
% (one-sample Kolmogorov-Smirnov tests, separately for each time series).
% When normal fits fail, the permutation results might be saved out into
% per-subject files depending on the value of "failedFitAction". The
% resulting files are named
% 'SUBJECTUMBER_FREQUENCYBAND_METHOD_failedFits.mat'.
%
% Assumes that the data is in EEGlab structures and that file naming
% follows the 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat')
% convention. Files are expected to be under a FREQUENCYBAND subfolder in
% "dirName" (e.g. [dirName, '/', FREQUENCYBAND, '/s05_alpha.mat']).
% 
% Further assumes that data is in not in analytic form, but 
% is real valued time series. For a surrogate estimation 
% function that expects analytic input (with separate angle and 
% envelope files) see surrEdgeEstimation.m.
%
% Also assumes that all relevant files are in the working directory.
%
% Optional input args are inferred from input arg types and values.
%
% Mandatory input:
% freq      - Frequency band (string) to work with. Needs to be the same 
%        as used in the file names (e.g. 'alpha' for files like 
%        's01_alpha.mat'). One of {'alpha', 'beta', 'gamma', 'delta',
%        'theta'}
% 
% Optional inputs:
% dirName   - Char array, path to folder containing the folder of 
%       frequency band data. E.g. subject-level data files if 
%       "freq" == 'alpha' are expected to be located under 
%       [dirName, '/', freq, '/']. Default is current working
%       directory (pwd).
% subjects  - Cell array, list of subjects whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
% method    - Either a char array, one of 
%       {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'}, or a cell array
%       of char arrays. Specifies one or more connectivity measures to
%       calculate on the surrogate data. Defaults to 'iplv'.
% dRate     - Decimation rate. We usually work with bandpass-filtered data
%       that nevertheless retains the original (1000 or 500 Hz) sampling
%       rate. For efficiency, we can decimate this data by providing a
%       dRate integer. E.g. dRate = 10 corresponds to one 10th of the data
%       remaining after decimation. Integer between 1-20, defaults to 1 
%       (no decimation). NOTE THAT THERE IS NO ADDITIONAL FILTERING STEP 
%       INCLUDED, BEWARE OF ALIASING EFFECTS
% surrNo    - Numeric value, number of surrogate data sets generated for
%       statistical testing of edge values. One of 100:100:20000, defaults 
%       to 10^3. 
% truncated - Either a char array, one of {'nontruncated', 'truncated'}, 
%       or a cell array of char arrays, one for each connectivity method 
%       specified in input arg "method". Determines if a standard or
%       truncated normal distribution should be fitted to surrogate edge 
%       data for the corresponding connectivity method. If 'truncated', 
%       a truncated normal is fitted (with bounds [0 1]), corresponding 
%       to the range of the phase-based connectivity measures. Defaults to 
%       'truncated' for phase-based connectivity measures and to 
%       'nontruncated' for envelope correlation. The latter is 
%       Fisher-Z transformed before fitting.
%       NOTE THAT THE RANGE [0 1] IS CURRENTLY HARDCODED!
% failedFitAction   - Char array, one of {'saveResults', 'nosave'}. Flag
%       defining how the function should behave upon encountering a failed 
%       (truncated) normal fit to the surrogate data. If set to 
%       'saveResults', per-epoch files are saved out, containing the 
%       permutation results of each failed fit. If 'nosave', 
%       fails are ignored. Files are named 
%       'SUBJECTUMBER_FREQUENCYBAND_METHOD_epochEPOCHNO_failedFits.mat' and
%       are saved out into a subfolder named SUBJECTNUMBER_failedFits in
%       the same folder where data files are located.
%       Defaults to 'saveResults'.
% 
% Output:
% The following variables are saved out for each subject.
% surrNormalMu      - Numeric array, sized 
%               (methods X ROIs X ROIs X epochs X conditions). Contains 
%               the estimated "mu" param for the normal / truncated normal 
%               distribution fitted to the surrogate connectivity data.
% surrNormalSigma   - Numeric array, sized 
%               (methods X ROIs X ROIs X epochs X conditions). Contains 
%               the estimated "sigma" param for the normal / truncated normal 
%               distribution fitted to the surrogate connectivity data.
% surrNormalH       - Numeric array, sized 
%               (methods X ROIs X ROIs X epochs X conditions). Elements 
%               are either 0 or 1. Contains the main result of the 
%               Kolmogorov-Smirnov goodness-of-fit test (kstest) testing 
%               if the surrogate connectivity data for each edge indeed 
%               corresponds to the fitted (truncated) normal distribution. 
%               1 = rejected at 0.05 level, 0 = null hypothesis cannot be 
%               rejected.
% surrNormalP       - Numeric array, sized 
%               (methods X ROIs X ROIs X epochs X conditions). Elements 
%               are in the range [0 1]. Contains the probability that the
%               surrogate connectivity data corresponds to the fitted 
%               (truncated) normal distribution. Outcome of the 
%               Kolmogorov-Smirnov goodness-of-fit test (kstest).
% method            - Cell array, contains the value(s) of input arg "method".
% surrNo            - Numeric value, value of input arg "surrNo".
% dRate             - Numeric value, value of input arg "dRate".
% subject           - Char array, content of input arg "subjects"
%               corresponding to the subject whose output file this var is
%               saved out to.
%
%
% NOTES: 
% (1) % The method for truncated normal is from:
% https://www.mathworks.com/matlabcentral/fileexchange/64040-fitting-a-truncated-normal-gaussian-distribution
% Alexey Ryabov
%
% TODO:
% - switch to EEGLAB filters
%


%% Input checks

% check number of args
if ~ismembertol(nargin, 1:8)
    error(['Function surrEdgeEstimationReal requires input arg "freq" ',...
        '(frequency band) while args "dirName", "subjects", "method", ',...
        '"dRate", "surrNo", "truncated" and "failedFitAction" are optional!']);
end
% check for mandatory argument
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end
% check optional arguments
% remember, "method" and "truncated" could be char arrays or cell arrays of
% char arrays as well
if ~isempty(varargin)
    for v = 1:length(varargin)    
        if iscell(varargin{v}) && ~all(ismember(varargin{v}, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'})) &&...
                ~all(ismember(varargin{v}, {'truncated', 'nontruncated'})) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif ischar(varargin{v}) && ~exist('dirName', 'var') && exist(varargin{v}, 'dir')
            dirName = varargin{v};
        elseif ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'})
            method = varargin{v};     
        elseif iscell(varargin{v}) && all(ismember(varargin{v}, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'})) && ~exist('method', 'var')
            method = varargin{v};            
        elseif isnumeric(varargin{v}) && ~exist('dRate', 'var') && ismember(varargin{v}, 1:20)
            dRate = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 100:100:20000)
            surrNo = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'truncated', 'nontruncated'}) && ~exist('truncated', 'var')
            truncated = varargin{v};  
        elseif iscell(varargin{v}) && all(ismember(varargin{v}, {'truncated', 'nontruncated'})) && ~exist('truncated', 'var')
            truncated = varargin{v};             
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'saveResults', 'nosave'}) && ~exist('failedFitAction', 'var')
            failedFitAction = varargin{v};
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "dirName", "subjects", "method", "dRate", ',...
                '"surrNo", "truncated" and "failedFitAction"!']);
        end
    end
end
% check if defaults are needed for input args
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
if ~exist('dRate', 'var')
    dRate = 1;
end
if ~exist('surrNo', 'var')
    surrNo = 10^3;
end
if ~exist('truncated', 'var')
    truncated = 'truncated';
end
if ~exist('failedFitAction', 'var')
    failedFitAction = 'saveResults';
end

% make sure that "method" is a cell, transform if char array - it is easier
% to treat it if the type is consistent
if ischar(method)
    method = {method};
end

% transform truncated to vector of logicals
% if there is only one value (char array), set to logical and expand it to
% match the number of measures specified in "method"
if ischar(truncated)
    if strcmp(truncated, 'nontruncated')
        truncated = false;
    elseif strcmp(truncated, 'truncated')
        truncated = true;
    end
    if numel(method) > 1
        truncated = repmat(truncated, [1, numel(method)]);
    end
% if "truncated" is cell, transform to boolean vector
elseif iscell(truncated)
    tmp = zeros(1, length(truncated));
    tmp(strcmp(truncated, 'truncated')) = 1;    
    truncated = logical(tmp);
end

% check if the length of "method" and "truncated" are the same, that is,
% there is one "truncated" flag for each connectivity "method"
if ~isequal(numel(method), numel(truncated))
    error('The length of input args "method" and "truncated" should be the same!');
end

% get char array versions of "method" to display for user and in file name
methodToDisplay = join(method, ', ');
methodToDisplay = methodToDisplay{:};
methodToFilename = join(method, '_');
methodToFilename = methodToFilename{:};

% issue warning if envelope correlation method is paired with truncated normals
if any(truncated & ismember(method, {'ampCorr', 'orthAmpCorr'}))
    warning(['Envelope correlation methods should be paired with a non-truncated normal fit. ',...
        'The function proceeds as you intend(?), but please double-check your use case!']);
end

% IMPORTANT!!! HARDCODED SAMPLE BOUNDS TO CONSIDER
sampleBounds = [1001 9000];

% user message
disp([char(10), char(10), 'IGNORE WARNINGS ABOUT TEMPORARY VARIABLES IN PARFOR LOOPS!', char(10), char(10)]);
disp([char(10), 'Starting surrEdgeEstimationReal_hyperscan4D function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure(s): ', methodToDisplay,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
    char(10), 'Truncated normal at [0 1]?: ', num2str(truncated), ...
    char(10), 'Action to take when fitting fails: ', failedFitAction,...
    char(10), 'Subjects: ']);
disp(subjects);
disp(['NOTE THAT THIS VERSION ONLY CONSIDERS SAMPLES BETWEEN ',... 
    num2str(sampleBounds(1)), ' AND ', num2str(sampleBounds(2)), '!!!']);
disp('NOTE THAT THIS VERSION EXPECTS 4D DATA!!!');


%% Basics: set params, functions, prepare for parfor loop

% get boolean flag from failedFitAction
if strcmp(failedFitAction, 'saveResults')
    failedFitSave = true;
end

% number  of connectivity methods
methodNo = length(method);

% number of subjects
subNo = length(subjects);

% standardize dirName format with regards to last char
if dirName(end) == '/'
    dirName = dirName(1:end-1);
end

% data folder
dataDir = [dirName, '/', freq, '/'];

% list all files we will use
dataFiles = cell(subNo, 1);
for s = 1:subNo
    dataFiles{s} = [dataDir, subjects{s}, '_', freq, '.mat'];
end

% load first data file, check dimensions
disp([char(10), 'Size check on the first data set...']);
subData = load(dataFiles{1});
subData = subData.EEG.data(:, sampleBounds(1):sampleBounds(2), :, :);
if numel(size(subData)) ~= 4
    error(['Expected 4D array for subject data, received ',... 
        num2str(numel(size(subData))), 'D data. Investigate!']);
else
    [roiNo, sampleNo, epochNo, condNo] = size(subData);
end

% if decimation is requested, update sample number accordingly
if dRate ~= 1
    subData = subData(:, 1:dRate:end, :, :);
    sampleNoOrig = sampleNo;
    sampleNo = size(subData, 2);
end

% Prepare a lowpass filter if envelope correlation is used
% Due to parfor, we need to define lpFilter even for cases when it is not
% used
if ismember(method, {'ampCorr', 'orthAmpCorr'})
    lpFilter = designfilt('lowpassiir',... 
            'PassbandFrequency', 10,...
            'StopbandFrequency', 13,... 
            'PassbandRipple', 0.1,... 
            'StopbandAttenuation', 50,... 
            'SampleRate', 1000/dRate,... 
            'MatchExactly',... 
            'passband');
else
    lpFilter = [];
end

% Define helper functions if truncated normal is to be fitted
% Due to parfor, we need to declare it even if it is not used (truncated is
% false)
%if truncated
% bounds - currently supported methods all have a range of [0 1]
x_min = 0;
x_max = 1;
% define heaviside as we would use it, with 1 at origin
heaviside_l = @(x) 1.0*(x>=0);
heaviside_r = @(x) 1.0*(x<=0);
% cdf of normal between bounds - needed for normalization
normcdf_lr =@(mu, sigma) (normcdf(x_max, mu, sigma) - normcdf(x_min, mu, sigma));
% truncated normal pdf
norm_trunc =@(x, mu, sigma) normpdf(x , mu, sigma)./normcdf_lr(mu, sigma) .* heaviside_l(x - x_min) .* heaviside_r(x - x_max);
%end

% setting options for fitting procedure
opts = statset('MaxIter', 100, 'MaxFunEval', 100, 'FunValCheck', 'off');

% user message
if dRate == 1
    disp([char(10), 'First data file had ', num2str(condNo), ' conditions, ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
        ' samples for each epoch. Assuming same channel and sample size for each data file.']);
else
    disp([char(10), 'First data file had ', num2str(condNo), ' conditions, ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(epochNo), ' epochs and ', num2str(sampleNoOrig), ... 
        ' samples for each epoch. Assuming same same channel and sample size for each data file.']);
    disp([char(10), 'After decimation (dRate = ', num2str(dRate), '), there ',...
        'are ', num2str(sampleNo), ' samples for each epoch.']);
end

% clock for timing the full function run
funcClock = tic;


%% Loop through subjects

parfor subIdx = 1:subNo
    
    % user message
    disp([char(10), 'Working on data from subject ', subjects{subIdx}]);
    
    % for timekeeping start a subject-level clock
    subClock = tic;
    
    % subject's EEG data
    subData = load(dataFiles{subIdx});
    subData = subData.EEG.data(:, sampleBounds(1):sampleBounds(2), :, :);
    if dRate ~= 1
        subData = subData(:, 1:dRate:end, :, :);
    end
    % check data size
    [subRoiNo, subSampleNo, subEpochNo, subCondNo] = size(subData);
    if ~isequal([subRoiNo, subSampleNo, subEpochNo, subCondNo], [roiNo, sampleNo, epochNo, condNo])
        error(['Data for subject ', subjects{subIdx}, ' has unexpected size, investigate!']);
    end 
    
    % preallocate result matrices
    surrNormalMu = nan(methodNo, roiNo, roiNo, subEpochNo, subCondNo);  % for storing "mu" of fitted normals
    surrNormalSigma = surrNormalMu;  % for storing "sigma" of fitted normals
    surrNormalP = surrNormalMu;  % for storing the p-value from kstest
    surrNormalH = surrNormalMu;  % for storing "h" (hypothesis test outcome) from kstest
    
    % set a flag to mark if a subject-level folder for failed fits data
    % exists - the flag starts as false and is set to true if a folder is
    % created later
    if failedFitSave
        failedFitDirFlag = false;
    end
    
    
    % loop through conditions
    for condIdx = 1:subCondNo
    
        % loop through epochs
        for epochIdx = 1:subEpochNo

            % preallocate cell array for capturing surrogate data for failed normal
            % fits
            if failedFitSave
                failedFits = cell(methodNo, roiNo, roiNo);
                failedFitsCounter = 0;
            end       

            % generate surrogate datasets and calculate
            % connectivity matrices for them

            % for envelope correlations pass the lowpass filter object
            if ismember(method, {'ampCorr', 'orthAmpCorr'})
                surrConnData = getSurrConn(subData(:, :, epochIdx, condIdx), surrNo, method, lpFilter);
            else
                surrConnData = getSurrConn(subData(:, :, epochIdx, condIdx), surrNo, method);
            end

            % fit a normal distribution to each group of edge values,
            % procees per connectivity method
            for methodIdx = 1:methodNo
                % get current method
                currentMethod = method{methodIdx};

                % go through edge pairings
                for roi1 = 1:roiNo
                    for roi2 = 1:roiNo

                        % only go through edges in the upper triangle,
                        % excluding the diagonal as well
                        if roi2 > roi1

                            % get data for current ROI-pairing
                            surrData = squeeze(surrConnData(methodIdx, :, roi1, roi2));

                            % for envelope correlations, first use Fisher-Z
                            % transform on values
                            if ismember(currentMethod, {'ampCorr', 'orthAmpCorr'})
                                surrData = atanh(surrData);
                            end

                            % fitting

                            kstestFlag = 0;  % flag for proceeding with Kolmogorov-Smirnov test

                            % if truncated normal is to be fitted
                            if truncated(methodIdx)
                                % use try - catch in case fitting errors out (happens when param would reach Inf / NaN value)
                                try
                                    % fit truncated normal
                                    phat = mle(surrData , 'pdf', norm_trunc, 'start', [mean(surrData), std(surrData)], 'Options', opts); 
                                    % save out main params from fitted normal
                                    surrNormalMu(methodIdx, roi1, roi2, epochIdx, condIdx) = phat(1);
                                    surrNormalSigma(methodIdx, roi1, roi2, epochIdx, condIdx) = phat(2);
                                    % create a probability distribution for the
                                    % truncated normal with fitted params
                                    pd = makedist('normal', 'mu', phat(1), 'sigma', phat(2));
                                    pdToTest = truncate(pd, x_min, x_max);      
                                    kstestFlag = 1;
                                catch ME
                                    disp(['Fitting with truncated normal dist failed at subject ',... 
                                        num2str(subIdx), ', epoch ', num2str(epochIdx),...
                                        ', rois ', num2str(roi1), ' and ', num2str(roi2)]);
                                    disp(ME.message);
                                    % store permutation results
                                    if strcmp(failedFitAction, 'saveResults')
                                        failedFits{methodIdx, roi1, roi2} = surrData;
                                        failedFitsCounter = failedFitsCounter + 1;
                                    end
                                end  % try
                            % if standard normal is to be fitted, not truncated    
                            elseif ~truncated(methodIdx)
                                % use try - catch in case fitting errors out (happens when param would reach Inf / NaN value)
                                try
                                    % fit simple normal, note that fitdist
                                    % requires column vector as data input
                                    pdToTest = fitdist(surrData', 'normal');  % output is a prob.NormalDistribution object
                                    % save out main params from fitted normal
                                    surrNormalMu(methodIdx, roi1, roi2, epochIdx, condIdx) = pdToTest.mu;
                                    surrNormalSigma(methodIdx, roi1, roi2, epochIdx, condIdx) = pdToTest.sigma;      
                                    kstestFlag = 1;
                                catch ME
                                    disp(['Fitting with normal dist failed at subject ',... 
                                        num2str(subIdx), ', epoch ', num2str(epochIdx),...
                                        ', rois ', num2str(roi1), ' and ', num2str(roi2)]);
                                    disp(ME.message);
                                    % store permutation results
                                    if strcmp(failedFitAction, 'saveResults')
                                        failedFits{methodIdx, roi1, roi2} = surrData;
                                        failedFitsCounter = failedFitsCounter + 1;
                                    end
                                end  % try
                            end  % if truncated(methodIdx)

                            % test the goodness-of-fit with single sample
                            % Kolmogorov-Smirnov
                            if kstestFlag  % only if fitting succeeded
                                [surrNormalH(methodIdx, roi1, roi2, epochIdx, condIdx), surrNormalP(methodIdx, roi1, roi2, epochIdx, condIdx)] = kstest(surrData, 'CDF', pdToTest);
                            end

                        end  % if roi2 > roi1
                    end  % for roi2 loop
                end  % for roi1 loop     
            end  % for methodIdx loop

            % save out data for failed fits if flag is set
            if failedFitSave && failedFitsCounter ~= 0
                % first check if there is already a folder created for holding
                % failed fits data - create folder if it has not been made yet
                if ~failedFitDirFlag
                    failedFitDir = mkdir([dirName, '/' , freq, '/', subjects{subIdx}, '_failedFits']);
                    failedFitDirFlag = true;  % adjust flag
                end
                % save failed fits data for current epoch into a .mat file,
                % use the "matfile" method that can be used in a parfor loop
                saveFailedFits = [failedFitDir, '/',... 
                    subjects{subIdx}, '_', freq, '_', methodToFilename,... 
                    '_epoch', num2str(epochIdx), '_cond', num2str(condIdx),... 
                    '_failedFits.mat'];
                saveFF = matfile(saveFailedFits);
                saveFF.failedFits = failedFits;
                saveFF.failedFitsCounter = failedFitsCounter;
                saveFF.epochIdx = epochIdx;
            end


        end  % for epochIdx loop
    
    end  % for condIdx loop    

    % get a subject-specific file name for saving results
    saveF = [dirName, '/' , freq, '/', subjects{subIdx}, '_', freq, '_surrEdgeEstReal_', methodToFilename, '.mat'];

    % save results using matfile - this is mostly allowed in a parfor loop,
    % unlike the save command
    saveM = matfile(saveF);
    saveM.surrNormalMu = surrNormalMu;
    saveM.surrNormalSigma = surrNormalSigma;
    saveM.surrNormalP = surrNormalP;
    saveM.surrNormalH = surrNormalH;
    saveM.method = method;
    saveM.surrNo = surrNo;
    saveM.dRate = dRate;
    saveM.subject = subjects{subIdx};
    saveM.truncated = truncated;
    % include the bounds if truncated normals were used
    if truncated
        saveM.x_min_max = [x_min, x_max];
    end    
    
    % report elapsed time
    disp([char(10), 'Finished with subject ', subjects{subIdx},...
        ', it took ', num2str(toc(subClock), 3), ' secs']);
        
end  % subIdx for loop
                    

%% Saving, cleaning up

disp([char(10), 'Done! Took ', num2str(toc(funcClock), 3),... 
    ' secs for the whole set']);

return

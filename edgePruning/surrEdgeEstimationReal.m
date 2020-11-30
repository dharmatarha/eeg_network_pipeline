function surrEdgeEstimationReal(freq, varargin)

%% Estimate distribution of edge weights (adj. matrix) for surrogate / null data
% To be used on real data, not analytic.
%
% USAGE: surrEdgeEstimationReal(freq,
%                               dirName = pwd, 
%                               subjects = {'s02', 's03', ...}, 
%                               method = 'iplv',
%                               dRate = 1,
%                               surrNo = 1000,
%                               truncated = 'yes')
%
% Fits edge weights for phase-scrambling-based surrogate data with a 
% (truncated) normal distribution, estimating the parameters (mean and std). 
% USES PARFOR!
%
% Connectivity is measured by calling functions (plv, iplv or pli) defined 
% outside this script. IMPORTANT: We only support phase-based
% connectivity-measures as of now.
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'SUBJECTNUMBER_FREQUENCYBAND_surrEdgeEstReal.mat'.
%
% Output files contain the parameters of the fitted (truncated) normal
% distributions and the results of goodness-of-fit tests
% (one-sample Kolmogorov-Smirnov tests, separately for each time series).
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
% dirName   - Directory path (string) pointing to 'angle' and 'envelope' 
%       folders. Also used for creating a 'connectivity' folder for saving 
%       out results. Default is current working directory (pwd).
% subjects  - Cell array, list of subjects whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray defined in
%       "restingStateSubjects.mat" (var "subjectsRS").
% method    - Connectivity measure compatible with phase data, one of 
%       {'plv', 'iplv', 'pli'}. Defaults to 'iplv'.
% dRate     - Decimation rate. We usually work with bandpass-filtered data
%       that nevertheless retains the original (1000 or 500 Hz) sampling
%       rate. For efficiency, we can decimate this data by providing a
%       dRate integer. E.g. dRate = 10 corresponds to one 10th of the data
%       remaining after decimation. Integer between 1-20, defaults to 1 
%       (no decimation). NOTE THAT THERE IS NO ADDITIONAL FILTERING STEP 
%       INCLUDED, BEWARE OF ALIASING EFFECTS
% surrNo    - Number of surrogate data sets generated for statistical
%       testing of edge values. Num value, one of 100:100:20000, defaults 
%       to 10^3. 
% truncated - Char array, one of {'no', 'yes'}. Determines if a standard or
%       truncated normal distribution should be fitted to surrogate edge 
%       data. If 'yes', that is, a truncated normal is fitted, we
%       truncate with bounds [0 1], corresponding to the range of the 
%       phase-based connectivity measures. Defaults to 'yes'.
% 
% Output:
% The following variables are saved out for each subject.
% surrNormalMu      - Numeric array, sized (node no. X node no. X layer no. X
%               stim no.). Contains the estimated "mu" param for the normal
%               / truncated normal distribution fitted to the surrogate
%               connectivity data.
% surrNormalSigma   - Numeric array, sized (node no. X node no. X layer no. X
%               stim no.). Contains the estimated "sigma" param for the normal
%               / truncated normal distribution fitted to the surrogate
%               connectivity data.
% surrNormalH       - Numeric array, sized (node no. X node no. X layer no. X
%               stim no.). Elements are either 0 or 1. Contains the main 
%               result of the Kolmogorov-Smirnov goodness-of-fit test 
%               (kstest) testing if the surrogate connectivity data for 
%               each edge indeed corresponds to the fitted (truncated) 
%               normal distribution. 1 = rejected at 0.05 level, 0 = null
%               hypothesis cannot be rejected.
% surrNormalP       - Numeric array, sized (node no. X node no. X layer no. X
%               stim no.). Elements are in the range [0 1]. Contains the 
%               probability that the surrogate connectivity data 
%               corresponds to the fitted (truncated) normal distribution. 
%               Outcome of the Kolmogorov-Smirnov goodness-of-fit test 
%               (kstest).
% method            - Char array, value of input arg "method".
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


%% Input checks

% check for mandatory argument
if ~ismembertol(nargin, 1:7)
    error(['Function surrEdgeEstimationReal requires input arg "freq" ',...
        '(frequency band) while args "dirName", "subjects", "method", ',...
        '"dRate", "surrNo" and "truncated" are optional!']);
end
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end
% check optional arguments
if ~isempty(varargin) 
    disp(varargin);
    for v = 1:length(varargin)    
        if iscell(varargin{v}) && ~exist('subjects', 'var')
            subjects = varargin{v};
        elseif ischar(varargin{v}) && ~exist('dirName', 'var') && exist(varargin{v}, 'dir')
            dirName = varargin{v};
        elseif ischar(varargin{v}) && ~exist('method', 'var') && ismember(varargin{v}, {'pli', 'iplv', 'plv'})
            method = varargin{v};     
        elseif isnumeric(varargin{v}) && ~exist('dRate', 'var') && ismember(varargin{v}, 1:20)
            dRate = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('surrNo', 'var') && ismember(varargin{v}, 100:100:20000)
            surrNo = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'yes', 'no'}) && ~exist('truncated', 'var')
            truncated = varargin{v};            
        else
            error(['There are either too many input args or they are not ',...
                'mapping nicely to "dirName", "subjects", "method", "dRate" and "surrNo"!']);
        end
    end
end
% check if defaults are needed for input args
if ~exist('dirName', 'var')
    dirName = pwd;
end
if ~exist('subjects', 'var')
    if exist('restingStateSubjects.mat' ,'file')
        tmp = load('restingStateSubjects.mat');
        subjects = tmp.subjectsRS;
    else
        error('Couldn''t find subject list in "restingStateSubjects.mat", no subject list to work on!');
    end
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
    truncated = 'yes';
end
% transform truncated to logical
if strcmp(truncated, 'no')
    truncated = false;
elseif strcmp(truncated, 'yes')
    truncated = true;
end

% user message
disp([char(10), 'Starting surrEdgeEstimationReal function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
    char(10), 'Truncated normal at [0 1]?: ', num2str(truncated), ...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% standardize dirName format with regards to last char
if dirName(end) == '/'
    dirName = dirName(1:end-1);
end

% data folder
dataDir = [dirName, '/', freq, '/'];

% list all angle and envelope files we will use
dataFiles = cell(subNo, 1);
for s = 1:subNo
    dataFiles{s} = [dataDir, subjects{s}, '_', freq, '.mat'];
end

% load first data file, check dimensions
disp([char(10), 'Size check on the first data set...']);
subData = load(dataFiles{1});
subData = subData.EEG.data;
if numel(size(subData)) ~= 3
    error(['Expected 3D array for subject data, received ',... 
        num2str(numel(size(subData))), 'D data. Investigate!']);
else
    [roiNo, sampleNo, epochNo] = size(subData);
end

% if decimation, update sample number accordingly
if dRate ~= 1
    subData = subData(:, 1:dRate:end, :, :);
    sampleNoOrig = sampleNo;
    sampleNo = size(subData, 2);
end

% define helper functions if truncated normal is to be fitted
if truncated
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
end

% user message
if dRate == 1
    disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
        ' samples for each epoch. Assuming same channel and sample size for each data file.']);
else
    disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
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
    
    % load angle and envelope data - for the first subject we have already 
    % loaded it in previous block, but due to how parfor works we need to
    % do so here again:
    
    % subject's EEG data
    subData = load(dataFiles{subIdx});
    subData = subData.EEG.data;
    if dRate ~= 1
        subData = subData(:, 1:dRate:end, :, :);
    end
    % check data size
    [sizeTmp1, sizeTmp2, ~] = size(subData);
    if ~isequal([sizeTmp1, sizeTmp2], [roiNo, sampleNo])
        error(['Data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
    end 
    
    % preallocate result matrices
    surrNormalMu = nan(roiNo, roiNo, epochNo);  % for storing "mu" of fitted normals
    surrNormalSigma = nan(roiNo, roiNo, epochNo);  % for storing "sigma" of fitted normals
    surrNormalP = nan(roiNo, roiNo, epochNo);  % for storing the p-value from kstest
    surrNormalH = nan(roiNo, roiNo, epochNo);  % for storing "h" (hypothesis test outcome) from kstest
    
   
    % loop through epochs
    for epochIdx = 1:epochNo

        % generate surrogate datasets and calculate
        % connectivity matrices for them
        surrConnData = getSurrConn(subData(:,:,epochIdx), surrNo, method);

        % fit a normal distribution to each group of edge values
        for roi1 = 1:roiNo
            for roi2 = 1:roiNo
                % only go through edges in the upper triangle,
                % excluding the diagonal as well
                if roi2 > roi1
                    tmp = squeeze(surrConnData(:, roi1, roi2));

                    % if truncated normal is to be fitted
                    if truncated
                        % fit truncated normal
                        phat = mle(tmp , 'pdf', norm_trunc, 'start', [mean(tmp), std(tmp)]);
                        % save out main params from fitted normal
                        surrNormalMu(roi1, roi2, epochIdx) = phat(1);
                        surrNormalSigma(roi1, roi2, epochIdx) = phat(2);
                        % create a probability distribution for the
                        % truncated normal with fitted params
                        pd = makedist('normal', 'mu', phat(1), 'sigma', phat(2));
                        pdToTest = truncate(pd, x_min, x_max);       
                    % if standard normal is to be fitted, not truncated    
                    else
                        % fitting
                        pdToTest = fitdist(tmp, 'normal');  % output is a prob.NormalDistribution object
                        % save out main params from fitted normal
                        surrNormalMu(roi1, roi2, epochIdx) = pd.mu;
                        surrNormalSigma(roi1, roi2, epochIdx) = pd.sigma;                     
                    end  % if truncated

                    % test the goodness-of-fit with single sample
                    % Kolmogorov-Smirnov
                    [surrNormalH(roi1, roi2, epochIdx), surrNormalP(roi1, roi2, epochIdx)] = kstest(tmp, 'CDF', pdToTest);

                end  % if roi2 > roi1
            end  % for roi2 loop
        end  % for roi1 loop          

    end  % for epochIdx  loop
    
    % get a subject-specific file name for saving results
    saveF = [dirName, '/' , subjects{subIdx}, '_', freq, '_surrEdgeEstReal.mat'];
    
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
    
    % report elapsed time
    disp([char(10), 'Finished with subject ', subjects{subIdx},...
        ', it took ', num2str(toc(subClock), 3), ' secs']);
        
end  % subIdx for loop
                    

%% Saving, cleaning up

disp([char(10), 'Done! Took ', num2str(toc(funcClock), 3),... 
    ' secs for the whole set']);

return
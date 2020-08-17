function surrEdgeEstimation(freq, varargin)

%% Estimate distribution of edge weights (adj. matrix) for surrogate / null data
%
% USAGE: surrEdgeEstimation(freq,
%                           dirName = pwd, 
%                           subjects = {'s02', 's03', ...}, 
%                           method = 'iplv',
%                           dRate = 1,
%                           surrNo = 1000,
%                           truncated = 'yes')
%
% Fits that edge weights for phase-scrambling-based surrogate data with a 
% normal distribution, estimating the parameters (mean and std). 
% USES PARFOR!
%
% Connectivity is measured by calling functions (plv, iplv or pli) defined 
% outside this script. IMPORTANT: We only support phase-based
% connectivity-measures as of now.
%
% Results are saved into the provided
% directory (or into current working directory), named
% 'SUBJECTNUMBER_FREQUENCYBAND_surrEdgeEstimate.mat'.
% Output files contain 
%
% Assumes that the data is in EEGlab structures and that file naming
% follows the 'SUBJECTNUMBER_FREQUENCYBAND.mat' (e.g. 's05_alpha.mat')
% convention
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
% subjects  - List of subjects (cell array) whose data we process. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
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
    error(['Function surrEdgeEstimation requires input arg "freq" ',...
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
            dirName = varargin{v};            
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
    subjects = {'s02','s03','s04','s05','s06','s07','s08','s09'...
     ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
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
    truncated = 'yes';
end
% transform truncated to logical
if strcmp(truncated, 'no')
    truncated = false;
elseif strcmp(truncated, 'yes')
    truncated = true;
end

% user message
disp([char(10), 'Starting surrEdgeEstimation function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
    char(10), 'Truncated normal at [0 Inf]?: ', num2str(truncated), ...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics

% number of subjects
subNo = length(subjects);

% angle and envelope folders
angleDir = [dirName, '/angle/', freq, '/'];
envDir = [dirName, '/envelope/', freq, '/'];

% list all angle and envelope files we will use
angleFiles = cell(subNo, 1);
envFiles = cell(subNo, 1);
for s = 1:subNo
    angleFiles{s} = [angleDir, subjects{s}, '_', freq, '.mat'];
    envFiles{s} = [envDir, subjects{s}, '_', freq, '.mat'];
end

% load first data files, both angle and envelope, check dimensions
disp([char(10), 'Sanity checks on the first angle and envelope data sets...']);
phaseData = load(angleFiles{1}, 'EEG');
phaseData = phaseData.EEG.data;
[roiNo, sampleNo, epochNo, stimNo] = size(phaseData);
% sanity check - is it phase data?
if any(any(any(any(phaseData>pi, 1), 2), 3), 4) || any(any(any(any(phaseData<-pi, 1), 2), 3), 4)
    error('There is data here outside the [-pi +pi] range, is it really phase data?');
end
% error if angle data has different size than envelope data
envData = load(envFiles{1}, 'EEG');
envData = envData.EEG.data;
if ~isequal([roiNo, sampleNo, epochNo, stimNo], size(envData))
    error('Angle and envelope data matrices have different dimensions, investigate!');
end

% decimation 
if dRate ~= 1
    phaseData = phaseData(:, 1:dRate:end, :, :);
    sampleNoOrig = sampleNo;
    sampleNo = size(phaseData, 2);
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
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNo), ... 
        ' samples for each epoch. Assuming same dimensions for each data file.']);
else
    disp([char(10), 'First data file had ', num2str(roiNo), ' channels/ROIs, ',... 
        num2str(stimNo), ' different stimuli (stories), ',...
        num2str(epochNo), ' epochs and ', num2str(sampleNoOrig), ... 
        ' samples for each epoch. Assuming same dimensions for each data file.']);
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
    
    % angle data
    phaseData = load(angleFiles{subIdx}, 'EEG');
    phaseData = phaseData.EEG.data;
    if dRate ~= 1
        phaseData = phaseData(:, 1:dRate:end, :, :);
    end
    % check data size
    if ~isequal(size(phaseData), [roiNo, sampleNo, epochNo, stimNo])
        error(['Angle data for subject ', subjects{subIndex}, ' has unexpected size, investigate!']);
    end
    % sanity check - is it really phase data, i.e. between -pi +pi
    if any(any(any(any(phaseData>pi, 1), 2), 3), 4) || any(any(any(any(phaseData<-pi, 1), 2), 3), 4)
        error('There is data here outside the [-pi +pi] range, is it really phase data?');
    end
    % envelope data
    envData = load(envFiles{subIdx}, 'EEG');
    envData = envData.EEG.data;    
    if dRate ~= 1
        envData = envData(:, 1:dRate:end, :, :);
    end
    % check data size
    if ~isequal(size(envData), [roiNo, sampleNo, epochNo, stimNo])
        error(['Envelope data for subject ', subjects{subIdx}, ' has unexpected size, investigate!']);
    end        
    
    % preallocate result matrices
    surrNormalMu = nan(roiNo, roiNo, epochNo, stimNo);  % for storing "mu" of fitted normals
    surrNormalSigma = nan(roiNo, roiNo, epochNo, stimNo);  % for storing "sigma" of fitted normals
    surrNormalP = nan(roiNo, roiNo, epochNo, stimNo);  % for storing the p-value from kstest
    surrNormalH = nan(roiNo, roiNo, epochNo, stimNo);  % for storing "h" (hypothesis test outcome) from kstest
    
    % loop through stimulus type / stories
    for stimIdx = 1:stimNo       
        % loop through epochs
        for epochIdx = 1:epochNo
            
            % get real-valued time series data from polar analytic
            [realData, ~] = pol2cart(squeeze(phaseData(:, :, epochIdx, stimIdx)), squeeze(envData(:, :, epochIdx, stimIdx))); 
                    
            % generate surrogate datasets and calculate
            % connectivity matrices for them
            surrConnData = getSurrConn(realData, surrNo, method);
            
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
                            surrNormalMu(roi1, roi2, epochIdx, stimIdx) = phat(1);
                            surrNormalSigma(roi1, roi2, epochIdx, stimIdx) = phat(2);
                            % create a probability distribution for the
                            % truncated normal with fitted params
                            pd = makedist('normal', 'mu', phat(1), 'sigma', phat(2));
                            pdToTest = truncate(pd, [x_min x_max]);       
                        % if standard normal is to be fitted, not truncated    
                        else
                            % fitting
                            pdToTest = fitdist(tmp, 'normal');  % output is a prob.NormalDistribution object
                            % save out main params from fitted normal
                            surrNormalMu(roi1, roi2, epochIdx, stimIdx) = pd.mu;
                            surrNormalSigma(roi1, roi2, epochIdx, stimIdx) = pd.sigma;                     
                        end  % if truncated
                        
                        % test the goodness-of-fit with single sample
                        % Kolmogorov-Smirnov
                        [surrNormalH(roi1, roi2, epochIdx, stimIdx), surrNormalP(roi1, roi2, epochIdx, stimIdx)] = kstest(tmp, 'CDF', pdToTest);
                        
                    end  % if roi2 > roi1
                end  % for roi2 loop
            end  % for roi1 loop          
            
        end  % for epochIdx  loop
    end  % for stimIdx loop
    
    % get a subject-specific file name for saving results
    saveF = [dirName, '/' , subjects{subIdx}, '_', freq, '_surrEdgeEstimate.mat'];
    
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
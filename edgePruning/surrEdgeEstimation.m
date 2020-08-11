function surrEdgeEstimation(freq, varargin)

%% Estimate distribution of edge weights (adj. matrix) for surrogate / null data
%
% USAGE: surrEdgeEstimation(freq, dirName = pwd, subjects = {'s02', 's03', ...}, method = 'iplv', dRate = 1, surrNo = 1000)
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
% 
%  


%% Input checks

% check for mandatory argument
if ~ismembertol(nargin, 1:6)
    error(['Function surrEdgeEstimation requires input arg "freq" ',...
        '(frequency band) while args "dirName", "subjects", "method", ',...
        '"dRate" and "surrNo" are optional!']);
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

% user message
disp([char(10), 'Starting surrEdgeEstimation function with following arguments: ',...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Working directory: ', dirName,...
    char(10), 'Connectivity measure: ', method,...
    char(10), 'Decimation rate: ', num2str(dRate),...
    char(10), 'No. of surrogate data sets: ', num2str(surrNo), ... 
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
                        % fitting
                        pd = fitdist(tmp, 'normal');  % output is a prob.NormalDistribution object
                        % test the goodness-of-fit with single sample
                        % Kolmogorov-Smirnov
                        [surrNormalH(roi1, roi2, epochIdx, stimIdx), surrNormalP(roi1, roi2, epochIdx, stimIdx)] = kstest(tmp, 'CDF', pd);
                        % save out main params from fitted normal
                        surrNormalMu(roi1, roi2, epochIdx, stimIdx) = pd.mu;
                        surrNormalSigma(roi1, roi2, epochIdx, stimIdx) = pd.sigma;
                    end  % if
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
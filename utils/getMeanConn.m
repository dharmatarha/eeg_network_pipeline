function meanConn = getMeanConn(freq, varargin)
%% Get group mean from connectivity data
%
% USAGE: meanConn = getMeanConn(freq, dirName=pwd, subjects={'s02', 's03', ...}, edgeThr=0, ) 
%
% Calculate mean connectivity matrices from a set of individual 
% connectivity results. 
% 
% Generally, we assume that individual connectivity results were
% generated with edgePruning.m. Accordingly, individual connectivity 
% results are assumed to be stored in files with the naming convention
% SUBJECTID_FREQUENCYBAND_edgePruningInfo.mat, e.g. 
% s02_alpha_edgePruningInfo.mat.
% Files are expected to be all located in "dirName", or in pwd if not
% provided. 
%
% Each individual connectivity result .mat contains a "realConn" var with
% the raw connectivity results, a "meanSurrConn" var with mean surrogate
% connectivity, a "prunedConn" var with pruned connectivity results and a
% "pValues" var with the significance values from real vs surrogate 
% comparisons. 
% Matrices "realConn" and "meanSurrConn" are averaged in a straightforward
% manner. "prunedConn" is expected to contain NaN values for 
% non-significant edges, the function first sets those values to zero. 
% Then the function averages across subjects by first selecting edges 
% present in more than "edgeThr" ratio of subjects (no selection takes place
% by default: edgeThr=0). 
%
% Optional arguments are sorted based on type.
%
% Mandatory input:
% freq      - Frequency band (string) to work with. Needs to be the same 
%        as used in the file names (e.g. 'alpha' for files like 
%        's01_alpha.mat'). One of {'alpha', 'beta', 'gamma', 'delta',
%        'theta'}
%
% Optional inputs:
% dirName   - Directory path (string) pointing to the data files. Also 
%       used for saving out group means. Default is current working 
%       directory (pwd).
% subjects  - List of subjects (cell array) whose data we average. 
%       Each cell contains a subject ID also used in the filenames 
%       (e.g. 's01' for files like 's01_alpha.mat'). If empty, 
%       we use the default cell arrray of:
%       subjects={'s02','s03','s04','s05','s06','s07','s08','s09'...
%           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
%           's21','s22','s23','s24','s25','s26','s27','s28'}
% edgeThr   - Numeric value determining the minimum ratio of individual
%       epochs and also of subjects that we expect averaged edges to be 
%       present in for mean pruned connectivity. Defaults to 0, must be
%       in range of 0:0.01:1.
%
% Output:
% meanConn  - Struct containing fields with averaged connectivity results
%       from "realConn", "meanSurrConn" and "prunedConn" matrices.
%


%% Input checks

% check number of arguments
if ~ismember(nargin, 1:4)
    error(['Function getMeanConn requires input "freq" and optional args ',...
        '"dirName", "subjects" and "edgeThr"!']);
end
% check mandatory arg
if ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" has an unexpected value!');
end
% sort varargin to optional args (if there is any)
if nargin > 1
    for v = 1:length(varargin)
        if ischar(varargin{v}) && exist(varargin{v}, 'dir')
            dirName = varargin{v};
        elseif iscell(varargin{v})
            subjects = varargin{v};
        elseif isnumeric(varargin{v}) && ismember(varargin{v}, 0:0.01:1)
            edgeThr = varargin{v};
        else
            error('There was an input that could not be mapped to any of "dirName", "subjects" or "edgeThr"!');
        end
    end
end
% define defaults if needed
if ~exist('dirName', 'var')
    dirName = pwd;
end
if ~exist('subjects', 'var')
    subjects = {'s02','s03','s04','s05','s06','s07','s08','s09'...
           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
           's21','s22','s23','s24','s25','s26','s27','s28'};
end
if ~exist('edgeThr', 'var')
    edgeThr = 0;
end

% user message
disp([char(10), 'Called getMeanConn function with input args: ', ...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Folder for files: ', dirName, ...
    char(10), 'Edge threshold ratio for pruned connectivity: ', num2str(edgeThr),...
    char(10), 'Subjects: ']);
disp(subjects);


%% Basics, settings

subNo = length(subjects);

% load first data, check dimensions
checkFile = [dirName, '/', subjects{1}, '_', freq, '_edgePruningInfo.mat'];
load(checkFile);
[roiNoOne, roiNoTwo, epochNo, stimNo] = size(realConn);
% first two dimensions need to match - square matric for connectivity data
if roiNoOne~=roiNoTwo
    error(['Var "realConn" from file ', checkFile,... 
        ' has different sizes for the first two dimensions, ',...
        'seemingly contains stg else than connectivity data!']);
end
roiNo = roiNoOne;
% all connectivity matrices need to have same sizes
if ~isequal(size(realConn), size(meanSurrConn)) || ~isequal(size(realConn), size(prunedConn))
    error('Sizes of connectivity result matrices "realConn", "prunedConn" and "meanSurrConn" are incostistent!');
end

% preallocate results struct
meanConn = struct;
meanConn.realConn = nan(roiNo, roiNo, epochNo, stimNo, subNo);
meanConn.realConnAvg = nan(roiNo, roiNo, epochNo, stimNo);
meanConn.meanSurrConn = nan(roiNo, roiNo, epochNo, stimNo, subNo);
meanConn.meanSurrConnAvg = nan(roiNo, roiNo, epochNo, stimNo);
meanConn.prunedConn = nan(roiNo, roiNo, epochNo, stimNo, subNo);
meanConn.prunedConnMask = nan(roiNo, roiNo, epochNo, stimNo);
meanConn.prunedConnAvg = nan(roiNo, roiNo, epochNo, stimNo);

% user messages
disp([char(10), 'First connectivity matrix had ', num2str(roiNo), ' channels/ROIs, ',... 
    num2str(stimNo), ' conditions (stimuli) and ',...
    num2str(epochNo), ' epochs per condition. ',...
    'Assuming same dimensions for each connectivity matrix.']);
disp([char(10), 'Decimation rate for first file was ', num2str(dRate),... 
    ', the connectivity measure was ', method,... 
    '. Assuming same params for each subject''s results.']);


%% Load all files, collect data

for subIdx = 1:subNo
    % load file into struct
    loadFile = [dirName, '/', subjects{subIdx}, '_', freq, '_edgePruningInfo.mat'];
    l = load(loadFile);
    % sanity checks
    if ~isequal(l.dRate, dRate)
        error(['File ', loadFile, ' had a different dRate than expected, investigate!']);
    end
    if ~isequal(l.method, method)
        error(['File ', loadFile, ' had a different connectivity method than expected, investigate!']);
    end        
    if ~isequal(size(l.realConn), [roiNo, roiNo, epochNo, stimNo]) || ~isequal(size(l.realConn), size(l.prunedConn)) || ~isequal(size(l.realConn), size(l.meanSurrConn))
        error(['File ', loadFile, ' had connectivity matrices with unexpected sizes, investigate!']);
    end
    % store full matrices
    meanConn.realConn(:, :, :, :, subIdx) = l.realConn;
    meanConn.prunedConn(:, :, :, :, subIdx) = l.prunedConn;
    meanConn.meanSurrConn(:, :, :, :, subIdx) = l.meanSurrConn;
end

% user message
disp([char(10), 'Loaded all data and collected individual matrices into struct "meanConn"']);


%% Averaging

% simple averages for "realConn" and "meanSurrNo"
meanConn.realConnAvg = mean(meanConn.realConn, 5);
meanConn.meanSurrConnAvg = mean(meanConn.meanSurrConn, 5);

% for "prunedConn" we create a mask first then select edges based on
% edgeThr
prunedConnMask = ~isnan(meanConn.prunedConn);  % mask of non-NaN edges
maskSum = sum(prunedConnMask, 5);  % sum across subjects
meanConn.prunedConnMask = maskSum > round(edgeThr*subNo);  % mask of edges present in more than edgeThr ratio of subjects
maskRepmat = repmat(meanConn.prunedConnMask, [1 1 1 1 subNo]);  % mask repeated for each subject, concatenated

% Edges not passing the threshold are set to NaN, similar to any edge that
% was not present in any individual data. Then all NaN values are set to 0
% for averaging - practically weighting edges with their prevalence across
% subjects
tmp = meanConn.prunedConn;  
tmp(~maskRepmat) = NaN;  % set edges not passing the threshold to NaN
tmp(isnan(tmp)) = 0;  % set all NaN values to 0
meanConn.prunedConnAvg = mean(tmp, 5);  % simple averaging

% set zero values back to NaN after averaging, to keep consistent with
% individual-level connectivity result files (stat scripts will except NaN
% values for non-sign. edges in pruned connectivity matrices)
meanConn.prunedConnAvg(meanConn.prunedConnAvg==0) = NaN; 

% user message
disp([char(10), 'Done, averaged all matrices']);


return










function mainPermComp(freq, method)

%% Helper / wrapper to calculate the permutation tests 
% on listener-listeren hyperscan data
%
% USAGE: mainPermComp(freq, method)
%
%
% Logic: 
% - For each frequency band, connectivity method and thresholding 
% (yes vs no) do:
%   - Group level:
%       - Permutation test comparing within- vs across-narrative similarity
%       - Permutation tests comparing within- vs across-narrative separately for
%       each narrative
%       - Permutation tests comparing within-narratives themselves with each
%       other
%   - Individual level:
%       - Permutation test comparing within- vs. across-narrative
%       similarity
% - Save out results, per freq and method


%% Input checks

if nargin ~= 2
    error('Requires input args "freq" and "method"!');
end
if ~ischar(freq) || ~ismember(freq, {'delta', 'theta', 'alpha', 'beta', 'gamma'})
    error('Input arg "freq" should be one of the main freq bands!');
end
if ~ischar(method) || ~ismember(method, {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "method" should be one of "plv", "iplv", "ampCorr" or "orthAmpCorr"!');
end


%% Basics, settings

% measure elapsed time for whole script
runClock = tic;

% base folder
% baseDir = '/media/adamb/bonczData/hyperscan/newSurrEdgeEstimates/';
baseDir = '/media/data_disk/adamb/speaker_listener_data/realdata/';

% % select connectivity data
% freq = 'alpha';
% method = 'plv';

% no. of subjects
subNo = 24;

% similarity metric
simMethod = 'corr';

% user message
disp([char(10), 'Requested comparison: ', ...
    char(10), 'Frequency band: ', freq, ...
    char(10), 'Connectivity measure: ', method, ...
    char(10), 'Similarity metric: ', simMethod]);


% save path
saveP = fullfile(baseDir, freq, ['mainCmpRes_', freq, '_', method, '.mat']);

% get file path
fileP = fullfile(baseDir, freq, ['group_surrResults_', freq, '_', method, '.mat']);
if ~exist(fileP, 'file')
    error(['Cannot find file at: ', fileP]);
end


%% Load data

data = load(fileP);


%% Group-level comparisons: Within- vs across-narrative similarity +
% within- vs across-narrative separately for each narrative +
% within-narratives themselves with each other

% results are collected into structs
groupRes_unthr = struct;
groupRes_thrSub = struct;
groupRes_thrGroup = struct;

% unthresholded
connData = permute(data.meanConn, [3 4 2 1]);
[groupRes_unthr.permRes,... 
    groupRes_unthr.withinCondPermRes,... 
    groupRes_unthr.connSim] = cmpFullConn(connData, simMethod, 10000, 'mean');

% thresholded - averaged after subject-specific thresholding
try
    connData = permute(data.meanMaskedConnPos+data.meanMaskedConnNeg, [3 4 2 1]);
    [groupRes_thrSub.permRes,... 
        groupRes_thrSub.withinCondPermRes,... 
        groupRes_thrSub.connSim] = cmpFullConn(connData, simMethod, 10000, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrSub.permRes = nan(1,1);
    groupRes_thrSub.withinCondPermRes = nan(1,1);
    groupRes_thrSub.connSim = nan(1,1);
end

% thresholded - thresholding on the group average itself
connData = permute(data.groupMaskedConnPos+data.groupMaskedConnNeg, [3 4 2 1]);
try
    [groupRes_thrGroup.permRes,... 
        groupRes_thrGroup.withinCondPermRes,... 
        groupRes_thrGroup.connSim] = cmpFullConn(connData, simMethod, 10000, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrGroup.permRes = nan(1,1);
    groupRes_thrGroup.withinCondPermRes = nan(1,1);
    groupRes_thrGroup.connSim = nan(1,1);
end
    
% display most interesting results
disp([char(10), char(10), 'Raw group average:']);
disp(groupRes_unthr.permRes(5));

if isstruct(groupRes_thrSub.permRes)
    disp([char(10), char(10), 'Group average of subject-thresholded data:']);
    disp(groupRes_thrSub.permRes(5));
end

if isstruct(groupRes_thrGroup.permRes)
    disp([char(10), char(10), 'Group-thresholded connectivity:']);
    disp(groupRes_thrGroup.permRes(5));
end


%% Subject-specific comparisons

% results are collected into structs
subRes_unthr = struct;
subRes_thr = struct;

for subIdx = 1:subNo
    
    % unthresholded
    connData = permute(squeeze(data.realConn(subIdx,:,:,:,:)), [3 4 2 1]);
    [subRes_unthr(subIdx).permRes,... 
        subRes_unthr(subIdx).withinCondPermRes,... 
        subRes_unthr(subIdx).connSim] = cmpFullConn(connData, simMethod, 10000, 'mean');

    % thresholded
    connData = permute(squeeze(data.maskedConnPos(subIdx,:,:,:,:)+data.maskedConnNeg(subIdx,:,:,:,:)), [3 4 2 1]);
    try
        [subRes_thr(subIdx).permRes,... 
            subRes_thr(subIdx).withinCondPermRes,... 
            subRes_thr(subIdx).connSim] = cmpFullConn(connData, simMethod, 10000, 'mean');
    catch ME
        disp(['Subject-level thresholded data calculation failed at subject ', num2str(subIdx)]);
        disp(ME)
        subRes_thr(subIdx).permRes = nan(1,1);
        subRes_thr(subIdx).withinCondPermRes = nan(1,1);
        subRes_thr(subIdx).connSim = nan(1,1);
    end
    
end


%% Save out results


save(saveP, 'groupRes_unthr', 'groupRes_thrSub', 'groupRes_thrGroup',...
    'subRes_unthr', 'subRes_thr',...
    'freq', 'method', 'simMethod'); 

% get elapsed time
et = round(toc(runClock), 3);
disp([char(10), 'Elapsed time for all comparisons: ', num2str(et), ' secs']);
disp('Done-done, bye-bye!');






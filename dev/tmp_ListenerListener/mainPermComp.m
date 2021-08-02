
function mainPermComp(freq, method)

%% Helper / wrapper to calculate the permutation tests 
% on listener-listeren hyperscan data
%
% USAGE: mainPermComp(freq, method)
%
%
% Logic: 
% 
% For given frequency band and connectivity method do:
% (1) Load surrogate stats, v1 ("group_surrResults_*.mat" file)
% (2) Get connectivity data with following thresholding options:
%       - Unthresholded
%       - Thresholded on subject level, separately for each epoch
%       - Thresholded on group level, separately for each epoch
% (3) Do connectivity similarity comparisons for each type
% (4) Repeat the above steps with surrogate stats v2 ("group_surrResultsv2_*.mat" file)
% (5) Repeat all of the above, but only take into account positively
% different edges (edges with connectivity > surrogate connectivity)
%
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
if ~ischar(method) || ~ismember(method, {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "method" should be one of "plv", "iplv", "ciplv", "ampCorr" or "orthAmpCorr"!');
end


%% Basics, settings

% measure elapsed time for whole script
runClock = tic;

% base folder
baseDir = '/media/adamb/bonczData/hyperscan/newSurrEdgeEstimates/';
% baseDir = '/media/data_disk/adamb/speaker_listener_data/realdata/';

% no. of subjects
if ismember(freq, {'theta', 'alpha', 'beta', 'gamma'})
    subNo = 26;
elseif strcmp(freq, 'delta')
    subNo = 25;
end

% similarity metric
simMethod = 'corr';

% permutations for random permutation test
permNo = 10000;

% user message
disp([char(10), 'Requested comparison: ', ...
    char(10), 'Frequency band: ', freq, ...
    char(10), 'Connectivity measure: ', method, ...
    char(10), 'Similarity metric: ', simMethod]);


% save path
saveP = fullfile(baseDir, freq, ['mainCmpRes_', freq, '_', method, '.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1 surrogate stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data

% get file path for surrogate stats
if strcmp(freq, 'gamma') && strcmp(method, 'plv') && contains(baseDir, 'bonczData')
    fileP = fullfile(baseDir, freq, ['group_surrResults_', freq, '_', method, '_True.mat']);
else
    fileP = fullfile(baseDir, freq, ['group_surrResults_', freq, '_', method, '.mat']);
end
if ~exist(fileP, 'file')
    error(['Cannot find file at: ', fileP]);
end

data = load(fileP);


%% Group-level comparisons: Within- vs across-narrative similarity +
% within- vs across-narrative separately for each narrative +
% within-narratives themselves with each other

% results are collected into structs
groupRes_unthr = struct;
groupRes_thrSub = struct;
groupRes_thrSub_pos = struct;
groupRes_thrGroup = struct;
groupRes_thrGroup_pos = struct;

% unthresholded
connData = permute(data.meanConn, [3 4 2 1]);
[groupRes_unthr.permRes,... 
    groupRes_unthr.withinCondPermRes,... 
    groupRes_unthr.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');

% thresholded - averaged after subject-specific thresholding, both negative
% and positive edge differences
try
    connData = permute(data.meanMaskedConnPos+data.meanMaskedConnNeg, [3 4 2 1]);
    [groupRes_thrSub.permRes,... 
        groupRes_thrSub.withinCondPermRes,... 
        groupRes_thrSub.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrSub.permRes = nan(1,1);
    groupRes_thrSub.withinCondPermRes = nan(1,1);
    groupRes_thrSub.connSim = nan(1,1);
end

% thresholded - averaged after subject-specific thresholding, only 
% positive edge differences
try
    connData = permute(data.meanMaskedConnPos, [3 4 2 1]);
    [groupRes_thrSub_pos.permRes,... 
        groupRes_thrSub_pos.withinCondPermRes,... 
        groupRes_thrSub_pos.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrSub_pos.permRes = nan(1,1);
    groupRes_thrSub_pos.withinCondPermRes = nan(1,1);
    groupRes_thrSub_pos.connSim = nan(1,1);
end

% thresholded - thresholding on the group average itself, both negative
% and positive edge differences
connData = permute(data.groupMaskedConnPos+data.groupMaskedConnNeg, [3 4 2 1]);
try
    [groupRes_thrGroup.permRes,... 
        groupRes_thrGroup.withinCondPermRes,... 
        groupRes_thrGroup.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrGroup.permRes = nan(1,1);
    groupRes_thrGroup.withinCondPermRes = nan(1,1);
    groupRes_thrGroup.connSim = nan(1,1);
end
    
% thresholded - thresholding on the group average itself, only 
% positive edge differences
connData = permute(data.groupMaskedConnPos, [3 4 2 1]);
try
    [groupRes_thrGroup_pos.permRes,... 
        groupRes_thrGroup_pos.withinCondPermRes,... 
        groupRes_thrGroup_pos.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_thrGroup_pos.permRes = nan(1,1);
    groupRes_thrGroup_pos.withinCondPermRes = nan(1,1);
    groupRes_thrGroup_pos.connSim = nan(1,1);
end


% display most interesting results
disp([char(10), char(10), 'VERSION 1 SURROGATE TESTS!!!']);

disp([char(10), 'Raw group average:']);
disp(groupRes_unthr.permRes(5));

disp([char(10), 'Group average of subject-thresholded data:']);
disp(groupRes_thrSub.permRes(5));

disp([char(10), 'Group average of subject-thresholded data, only positive edges:']);
disp(groupRes_thrSub_pos.permRes(5));

disp([char(10), 'Group-thresholded connectivity:']);
disp(groupRes_thrGroup.permRes(5));

disp([char(10), 'Group-thresholded connectivity, only positive edges:']);
disp(groupRes_thrGroup_pos.permRes(5));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V2 surrogate stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data

% get file path for surrogate stats
fileP = fullfile(baseDir, freq, ['group_surrResultsv2_', freq, '_', method, '.mat']);
if ~exist(fileP, 'file')
    error(['Cannot find file at: ', fileP]);
end

data = load(fileP);


%% Group-level comparisons: Within- vs across-narrative similarity +
% within- vs across-narrative separately for each narrative +
% within-narratives themselves with each other

% results are collected into structs
groupRes_v2_thrSub = struct;
groupRes_v2_thrSub_pos = struct;
groupRes_v2_thrGroup = struct;
groupRes_v2_thrGroup_pos = struct;

% thresholded - averaged after subject-specific thresholding, both negative
% and positive edge differences
try
    connData = permute(data.meanMaskedConnPos+data.meanMaskedConnNeg, [3 4 2 1]);
    [groupRes_v2_thrSub.permRes,... 
        groupRes_v2_thrSub.withinCondPermRes,... 
        groupRes_v2_thrSub.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_v2_thrSub.permRes = nan(1,1);
    groupRes_v2_thrSub.withinCondPermRes = nan(1,1);
    groupRes_v2_thrSub.connSim = nan(1,1);
end

% thresholded - averaged after subject-specific thresholding, only 
% positive edge differences
try
    connData = permute(data.meanMaskedConnPos, [3 4 2 1]);
    [groupRes_v2_thrSub_pos.permRes,... 
        groupRes_v2_thrSub_pos.withinCondPermRes,... 
        groupRes_v2_thrSub_pos.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_v2_thrSub_pos.permRes = nan(1,1);
    groupRes_v2_thrSub_pos.withinCondPermRes = nan(1,1);
    groupRes_v2_thrSub_pos.connSim = nan(1,1);
end

% thresholded - thresholding on the group average itself, both negative
% and positive edge differences
connData = permute(data.groupMaskedConnPos+data.groupMaskedConnNeg, [3 4 2 1]);
try
    [groupRes_v2_thrGroup.permRes,... 
        groupRes_v2_thrGroup.withinCondPermRes,... 
        groupRes_v2_thrGroup.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_v2_thrGroup.permRes = nan(1,1);
    groupRes_v2_thrGroup.withinCondPermRes = nan(1,1);
    groupRes_v2_thrGroup.connSim = nan(1,1);
end
    
% thresholded - thresholding on the group average itself, only 
% positive edge differences
connData = permute(data.groupMaskedConnPos, [3 4 2 1]);
try
    [groupRes_v2_thrGroup_pos.permRes,... 
        groupRes_v2_thrGroup_pos.withinCondPermRes,... 
        groupRes_v2_thrGroup_pos.connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
catch ME
    disp('Group-thresholded data calculation failed!');
    disp(ME)
    groupRes_v2_thrGroup_pos.permRes = nan(1,1);
    groupRes_v2_thrGroup_pos.withinCondPermRes = nan(1,1);
    groupRes_v2_thrGroup_pos.connSim = nan(1,1);
end


% display most interesting results

disp([char(10), char(10), 'VERSION 2 SURROGATE TESTS!!!']);

disp([char(10), 'Group average of subject-thresholded data:']);
disp(groupRes_v2_thrSub.permRes(5));

disp([char(10), 'Group average of subject-thresholded data, only positive edges:']);
disp(groupRes_v2_thrSub_pos.permRes(5));

disp([char(10), 'Group-thresholded connectivity:']);
disp(groupRes_v2_thrGroup.permRes(5));

disp([char(10), 'Group-thresholded connectivity, only positive edges:']);
disp(groupRes_v2_thrGroup_pos.permRes(5));


% 
% %% Subject-specific comparisons
% 
% % results are collected into structs
% subRes_unthr = struct;
% subRes_thr = struct;
% 
% for subIdx = 1:subNo
%     
%     % unthresholded
%     connData = permute(squeeze(data.realConn(subIdx,:,:,:,:)), [3 4 2 1]);
%     [subRes_unthr(subIdx).permRes,... 
%         subRes_unthr(subIdx).withinCondPermRes,... 
%         subRes_unthr(subIdx).connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
% 
%     % thresholded
%     connData = permute(squeeze(data.maskedConnPos(subIdx,:,:,:,:)+data.maskedConnNeg(subIdx,:,:,:,:)), [3 4 2 1]);
%     try
%         [subRes_thr(subIdx).permRes,... 
%             subRes_thr(subIdx).withinCondPermRes,... 
%             subRes_thr(subIdx).connSim] = cmpFullConn(connData, simMethod, permNo, 'mean');
%     catch ME
%         disp(['Subject-level thresholded data calculation failed at subject ', num2str(subIdx)]);
%         disp(ME)
%         subRes_thr(subIdx).permRes = nan(1,1);
%         subRes_thr(subIdx).withinCondPermRes = nan(1,1);
%         subRes_thr(subIdx).connSim = nan(1,1);
%     end
%     
% end


%% Save out results


save(saveP, 'groupRes_unthr', 'groupRes_thrSub', 'groupRes_thrSub_pos', ...
    'groupRes_thrGroup', 'groupRes_thrGroup_pos', ...
    'groupRes_v2_thrSub', 'groupRes_v2_thrSub_pos', ...
    'groupRes_v2_thrGroup', 'groupRes_v2_thrGroup_pos', ...
    'freq', 'method', 'simMethod'); 

% get elapsed time
et = round(toc(runClock), 3);
disp([char(10), 'Elapsed time for all comparisons: ', num2str(et), ' secs']);
disp('Done-done, bye-bye!');






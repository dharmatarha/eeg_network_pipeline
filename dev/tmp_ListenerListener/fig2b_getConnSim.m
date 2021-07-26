%% Script for Figure 2B
% The script generates the group-averaged connectivity similarity matrix 
% depicted in fig2B
%
%

%% Base params

method = 'plv'; 
freq = 'alpha';
simMethod = 'corr';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['mainCmpRes_', freq, '_', method, '.mat']);
res = load(fileP);

% sanity check, correct similarity method?
if ~strcmp(res.simMethod, simMethod)
    warning(['WRONG SIMILARITY METHOD? LOADED RESULTS FILE WAS OBTAINED WITH ', upper(res.simMethod), '!']);
end


%% Get connectivity similarity

connSim = res.(['groupRes_', thr]).connSim;

        
%% Plot as heatmap

% punch the diagonal out
connSim(1:size(connSim, 1)+1:end) = nan;

% plotting constants
gcfMainPos = [0, 0, 0.58, 96];
gcfBackground = [1 1 1];
missingEdgeLabel = 'NA';
missingDataColor = [0 0 0];
fontSize = 14;

% heatmap 
h = heatmap(connSim, 'ColorMap', jet,... 
    'MissingDataColor', missingDataColor,... 
    'MissingDataLabel', missingEdgeLabel,... 
    'FontSize', fontSize);

% background to white
set(gcf,'Color', gcfBackground);
% set size
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);

% remove tick labels - release specific method accessing private properties:
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);

hs.XAxis.TickValues = categorical([41, 81, 121]);
hs.YAxis.TickValues = categorical([41, 81, 121]);
% hs.XAxis.TickValues = labelsRo;
% hs.YAxis.TickValues = labelsRo;


%% Save out, quit

saveas(gcf, fullfile(baseDir, ['connSim_', thr, '_', freq, '_', method, '.svg']));
saveas(gcf, fullfile(baseDir, ['connSim_', thr, '_', freq, '_', method, '.png']));





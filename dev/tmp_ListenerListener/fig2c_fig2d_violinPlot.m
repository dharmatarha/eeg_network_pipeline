%% Script for Figure 2C and 2D
% The script generates the violin plots depicted in fig2C and fig2d.
%
%

%% Base params

method = 'ciplv'; 
freq = 'alpha';
metric = 'corr';  % similarity metric

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['group_surrResultsv2_', freq, '_', method, '.mat']);
res = load(fileP);


%% Get mean connectivity, depending on "thr"
% Only versions where the positively larger connections are considered

switch thr
    case 'unthr'
        connData = res.meanConn;
    case 'thrSub'
        % get the data from the averaged, subject-level thresholded matrices, only positive diffs
        connData = res.meanMaskedConnPos;   
    case 'thrGroup'
        connData = res.groupMaskedConnPos;
end  % switch


% reorder so that first two dimensions yield a connectivity matrix for
% given epoch and condition
connData = permute(connData, [3 4 2 1]);
% sanity check
if ~isequal(size(connData), [62, 62, 40, 4])
    error('Wrong connData var size!');
end


%% Get connectivity similarity

[roiNo, ~, epochNo, condNo] = size(connData);

% reshape data so that epochs across conditions/stimuli come after each
% other
connData = reshape(connData, [roiNo, roiNo, epochNo*condNo]);

% calculation depends on metric type

% simple correlation
if strcmp(metric, 'corr')
    % extract upper triangles above the main diagonal
    dataLin = linearizeTrius(connData, 1);  
    % get all column-pairwise correlations
    connSim = corr(dataLin, dataLin);
    % keep the upper triangle, inlcuding main diagonal, set the rest to NaN
    connSim(tril(true(epochNo*condNo), -1)) = nan;
    
% frobenius norm of difference matrix ('eucl') or DeltaCon
elseif ismember(metric, {'eucl', 'deltaCon'})
    % preallocate results matrix
    connSim = nan(epochNo*condNo);
    % loops through epochs
    for epochOne = 1:epochNo*condNo
        % create symmetric adjacency matrices with zeros at diagonal
        epochOneData = triu(connData(:, :, epochOne), 1) + triu(connData(:, :, epochOne), 1)';
        for epochTwo = 1:epochNo*condNo
            if epochOne <= epochTwo  % only for upper triangle of pairings
                
                % create symmetric adjacency matrices with zeros at diagonal
                epochTwoData = triu(connData(:, :, epochTwo), 1) + triu(connData(:, :, epochTwo), 1)';
                
                % calculation depending on metric
                switch metric
                    case 'eucl'
                        connSim(epochOne, epochTwo) = norm(epochOneData-epochTwoData, 'fro');
                    case 'deltaCon'
                        connSim(epochOne, epochTwo) = deltaCon(epochOneData, epochTwoData, false);  % verbosity of deltaCon is set to false
                end
                
            end  % if
        end  % for groupEpoch
    end  % for subEpoch      
    
end
        

%% Get within- vs across-condition similarities

% first mirror the upper triangular part of connSim matrix to the lower
% half
connSim = triu(connSim, 1) + triu(connSim, 1)'; 

% preallocate vars holding within- and across condition/stimulus similarities
withinCondSim = nan((epochNo^2-epochNo)/2, condNo);
acrossCondSim = nan(epochNo^2*(condNo-1), condNo);

for condIdx = 1:condNo
    
    % within-cond epoch pairing similarities for given condition
    tmp = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, (condIdx-1)*epochNo+1:condIdx*epochNo); 
    withinCondSim(:, condIdx) = tmp(triu(true(epochNo), 1));
    
    % get number of NaN values among similarity values - occurs e.g. from fully
    % zero connectivity matrices compared with "corr"
    withinCondNan = sum(isnan(withinCondSim(:, condIdx)));
    
    % across-condition pairings for epochs in given condition
    tmpParts = cell(2, 1);
    if condIdx ~= 1
        tmpParts{1} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, 1:(condIdx-1)*epochNo);
    end
    if condIdx ~= condNo
        tmpParts{2} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, condIdx*epochNo+1:end);
    end
    % concatenate the parts
    if isempty(tmpParts{1})
        tmpAll = tmpParts{2};
    elseif isempty(tmpParts{2})
        tmpAll = tmpParts{1};
    else
        tmpAll = [tmpParts{1}, tmpParts{2}];
    end
    % linearize into a vector
    acrossCondSim(:, condIdx) = reshape(tmpAll', [epochNo^2*(condNo-1), 1]);
    
    % get number of NaN values among similarity values - occurs e.g. from fully
    % zero connectivity matrices compared with "corr"
    acrossCondNan = sum(isnan(acrossCondSim(:, condIdx)));    
    
end


% %% Violin plot 2c - the four narratives in one plot
% 
% % data into cell array
% data = cell(1, 5); 
% for i = 1:4
%     data{i} = withinCondSim(:, i);
% end
% data{5} = acrossCondSim(:);
% % sanity check for NaN values
% for i=1:5
%     flag = any(isnan(data{i})); 
%     if flag
%         disp(['NaN value(s) in cell ', num2str(i), '!']);
%     end
% end
% 
% % plotting constants
% gcfMainPos = [0, 0, 0.40, 0.60];
% gcfBackground = [1 1 1];
% fontSize = 14;
% fontWeight = 'Bold';
% 
% % colors for violins
% facecolor = [0.8500, 0.3250, 0.0980;
%     0.8500, 0.3250, 0.0980;
%     0.8500, 0.3250, 0.0980;
%     0.8500, 0.3250, 0.0980];
% 
% % labels for violins
% labels = {'Context 1', 'Context 2', 'Context 3', 'Context 4'};
% 
% [h, l, mx, med, bw] = violin(data(1:4), 'medc', [], 'mc', 'k', 'xlabel', labels, ... 
%     'facealpha', 0.6, 'facecolor', facecolor, 'edgecolor', []);
% ylabel({'FC matrix similarity', 'across epochs (correlation)'});
% ylim([0.4, 0.8]);
% 
% % background to white
% set(gcf,'Color', gcfBackground);
% % set size
% set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
% 
% % Font size, weight
% set(gca, 'FontSize', fontSize);
% set(gca, 'FontWeight', fontWeight);
% 
% % Misc
% set(gca, 'YGrid', 'on');
% 
% saveas(gcf, fullfile(baseDir, ['narrativesViolins_', thr, '_', freq, '_', method, '.svg']));
% saveas(gcf, fullfile(baseDir, ['narrativesViolins_', thr, '_', freq, '_', method, '.png']));
% 
% 
% %% Violin plot 2d - the within- vs across-narratives data in one plot
% 
% close(gcf);
% 
% % data into cell array
% data = cell(1, 5); 
% for i = 1:4
%     data{i} = withinCondSim(:, i);
% end
% data{5} = acrossCondSim(:);
% % sanity check for NaN values
% for i=1:5
%     flag = any(isnan(data{i})); 
%     if flag
%         disp(['NaN value(s) in cell ', num2str(i), '!']);
%     end
% end
% 
% % data together for within-narratives
% dataWithin = cat(1, data{1}, data{2}, data{3}, data{4});
% dataAcross = data{5};
% 
% % plotting constants
% gcfMainPos = [0, 0, 0.30, 0.60];
% gcfBackground = [1 1 1];
% fontSize = 14;
% fontWeight = 'Bold';
% 
% % colors for violins
% facecolor = [0.8500, 0.3250, 0.0980;
%     0, 0.4470, 0.7410];
% 
% % labels for violins
% labels = {'Within-narratives', 'Across-narratives'};
% 
% [h, l, mx, med, bw] = violin({dataWithin, dataAcross}, 'medc', [], 'mc', 'k', 'xlabel', labels, ... 
%     'facealpha', 0.6, 'facecolor', facecolor, 'edgecolor', []);
% l.Location = 'north';
% ylabel({'Connectivity similarity', 'across epochs (correlation)'});
% ylim([0.4, 0.8]);
% 
% % background to white
% set(gcf,'Color', gcfBackground);
% % set size
% set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
% 
% % Font size, weight
% set(gca, 'FontSize', fontSize);
% set(gca, 'FontWeight', fontWeight);
% 
% % Misc
% set(gca, 'YGrid', 'on');
% 
% 
% saveas(gcf, fullfile(baseDir, ['withinVSacrossViolins_', thr, '_', freq, '_', method, '.svg']));
% saveas(gcf, fullfile(baseDir, ['withinVSacrossViolins_', thr, '_', freq, '_', method, '.png']));


%% Violin plot - fig 2c and 2d in one

% data into cell array
data = cell(1, 5); 
for i = 1:4
    data{i} = withinCondSim(:, i);
end
data{5} = acrossCondSim(:);
% sanity check for NaN values
for i=1:5
    flag = any(isnan(data{i})); 
    if flag
        disp(['NaN value(s) in cell ', num2str(i), '!']);
    end
end

% plotting constants
gcfMainPos = [0, 0, 0.58, 96];
gcfBackground = [1 1 1];
fontSize = 14;
fontWeight = 'Bold';

% colors for violins
facecolor = [0.8500, 0.3250, 0.0980;
    0.8500, 0.3250, 0.0980;
    0.8500, 0.3250, 0.0980;
    0.8500, 0.3250, 0.0980;
    0, 0.4470, 0.7410];

% labels for violins
labels = {'Context 1', 'Context 2', 'Context 3', 'Context 4', 'Across-context'};

[h, l, mx, med, bw] = violin(data, 'medc', [], 'mc', 'k', 'xlabel', labels, ... 
    'facealpha', 0.6, 'facecolor', facecolor, 'edgecolor', []);
% ylabel({'FC matrix similarity', 'across epochs (correlation)'});
ylabel('FC matrix similarity');
ylim([0.4, 0.8]);

% background to white
set(gcf,'Color', gcfBackground);
% set size
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);

% Font size, weight
set(gca, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight);

% Misc
set(gca, 'YGrid', 'on');
















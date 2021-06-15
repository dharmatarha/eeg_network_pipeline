
% % load zrand similarity results
% modZrandF = '/home/adamb/eeg_network_pipeline/mod_zrand_varGroupSize_alpha_orthAmpCorr.mat';
% load(modZrandF);
% 
% % similarity values in cell array
% y = cell(1, 11);
% y{1} = zrandSubjectsNorm;
% for i=2:11
%     y{i} = zrandResNorm(:, i-1);
% end
% 
% [h,L,MX,MED] = violin(y, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0.6350, 0.0780, 0.1840], 'facealpha', 1);
% 
% 
% xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
% xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
% xlabel('Number of subjects in groups');
% yticks([-0.2 0 0.2 0.4 0.6 0.8 1]);
% yticklabels({'-0.2', '0', '0.2', '0.4','0.6', '0.8', '1'});
% ylabel(['Between-group partition similarity', char(10), '(z-rand score normalized to [0 1])']);
% ylim([-0.2 1])
% hold on;
% connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
% connectingLineYvalues = MED;
% plot(connectingLineXvalues, connectingLineYvalues, 'k', 'LineWidth', 1);
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'w');
% grid on;
% 
% L.String = L.String(1);
% L.Location = 'southeast';



% base folder to work in
baseDir = '/media/adamb/bonczData/EEG_resting_state/';
% method
method = 'orthAmpCorr';

% load zrand similarity values for the subject-to-subject comparisons
subCompF = [baseDir, 'fullConnMod_alpha_', method, '_1.05.mat'];
tmp = load(subCompF);
zrandRes1 = tmp.zrandRes(:);

% load zrand similarity results for different group sizes
groupCompF = [baseDir, 'mod_zrand_varGroupSize_alpha_', method, '.mat'];
tmp = load(groupCompF);
zrandRes2 = tmp.zrandRes;

% normalize
maxValue = 43.5;
zrandSubjectsNorm = zrandRes1./maxValue;
zrandResNorm = zrandRes2./maxValue;

% similarity values in cell array
y = cell(1, 11);
y{1} = zrandSubjectsNorm;
for i=2:11
    y{i} = zrandResNorm(:, i-1);
end

% colors:
% plv: [0, 0.4470, 0.7410]
% iplv: [0.9290, 0.6940, 0.1250]
% ampCorr: [0.4660, 0.6740, 0.1880]
% orthAmpCorr: [0.6350, 0.0780, 0.1840]
switch method
    case 'plv'
        violinColor = [0, 0.4470, 0.7410];
        titleText = 'PLV';
    case 'iplv'
        violinColor = [0.9290, 0.6940, 0.1250];
        titleText = 'iPLV';
    case 'ampCorr'
        violinColor = [0.4660, 0.6740, 0.1880];
        titleText = 'ampCorr';
    case 'orthAmpCorr'
        violinColor = [0.6350, 0.0780, 0.1840];
        titleText = 'orthAmpCorr';
end

% call violin plotter
[h,L,MX,MED] = violin(y, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', violinColor, 'facealpha', 1);


xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
xlabel('Number of subjects in groups');
yticks([-0.2 0 0.2 0.4 0.6 0.8 1]);
yticklabels({'-0.2', '0', '0.2', '0.4','0.6', '0.8', '1'});
ylabel(['Between-group partition similarity', char(10), '(z-rand score normalized to [-1 1])']);
ylim([-0.2 1])
title(titleText);
hold on;
connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
connectingLineYvalues = MED;
plot(connectingLineXvalues, connectingLineYvalues, 'k', 'LineWidth', 1);
set(gca, 'FontSize', 11);
set(gcf, 'Color', 'w');
grid on;

L.String = L.String(1);
L.Location = 'southeast';

% load zrand similarity results
modZrandF = '/home/adamb/eeg_network_pipeline/mod_zrand_varGroupSize_alpha_orthAmpCorr.mat';
load(modZrandF);

% similarity values in cell array
y = cell(1, 11);
y{1} = zrandSubjectsNorm;
for i=2:11
    y{i} = zrandResNorm(:, i-1);
end

[h,L,MX,MED] = violin(y, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0.6350, 0.0780, 0.1840], 'facealpha', 1);


xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
xlabel('Number of subjects in groups');
yticks([-0.2 0 0.2 0.4 0.6 0.8 1]);
yticklabels({'-0.2', '0', '0.2', '0.4','0.6', '0.8', '1'});
ylabel(['Between-group partition similarity', char(10), '(z-rand score normalized to [0 1])']);
ylim([-0.2 1])
hold on;
connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
connectingLineYvalues = MED;
plot(connectingLineXvalues, connectingLineYvalues, 'k', 'LineWidth', 1);
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');
grid on;

L.String = L.String(1);
L.Location = 'southeast';

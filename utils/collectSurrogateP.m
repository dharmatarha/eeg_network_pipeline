

% load and aggregate surrogate edge distribution data files for given
% frequency
freq = 'alpha';
surrFolder = '/media/adamb/bonczData/hyperscan/surrEdgeEstimates/truncatedNormal/';
subIndices = [2:9, 11:28];

% preallocate
surrP = zeros(62*61/2*160, length(subIndices));  % variable holding all KS-test p value outcomes across all subjects
% counter for subjects
subCounter = 0;

for subIdx = subIndices

    % adjust counter
    subCounter = subCounter+1;
    
    % get file name
    if subIdx < 10
        surrSubFile = [surrFolder, freq, '/s0', num2str(subIdx),'_alpha_surrEdgeEstimate.mat'];
    else
        surrSubFile = [surrFolder, freq, '/s', num2str(subIdx),'_alpha_surrEdgeEstimate.mat'];
    end
    % load data
    subData = load(surrSubFile);
    % get p value for the fit of the truncated normal on estimates
    tmp = subData.surrNormalP;
    tmp(isnan(tmp)) = [];
    % store in aggregating var
    surrP(:, subCounter) = tmp;
    
end

% get ratio of tests rejected at different levels
totalTestNo = numel(surrP);
noRej1 = sum(sum(surrP<0.1));
noRej05 = sum(sum(surrP<0.05));
noRej01 = sum(sum(surrP<0.01));
noRej001 = sum(sum(surrP<0.001));
disp([newline, 'Tests rejected at p < 0.1: ', num2str(noRej1), ', a ratio of ', num2str(noRej1/totalTestNo), newline]);
disp([newline, 'Tests rejected at p < 0.05: ', num2str(noRej05), ', a ratio of ', num2str(noRej05/totalTestNo), newline]);
disp([newline, 'Tests rejected at p < 0.01: ', num2str(noRej01), ', a ratio of ', num2str(noRej01/totalTestNo), newline]);
disp([newline, 'Tests rejected at p < 0.001: ', num2str(noRej001), ', a ratio of ', num2str(noRej001/totalTestNo), newline]);
    











function modularityRealVersusNull(realResFile, nullResFile, xlim, ylim)
%% Generating heatmaps from corresponding outputs of multiCommDetectWrapper 
% and multiCommDetectWrapper runs.
%
% Inputs:
% realResFile   - Char array, path to file storing multiCommDetectWrapper
%               outputs
% nullResFile   - Char array, path to file storing multiCommDetectNullWrapper
%               outputs
% xlim          - Numeric vector, with numel(xlim)=2. X axis limits for
%               displaying results.
% ylim          - Numeric vector, with numel(ylim)=2. Y axis limits for
%               displaying results.
%
% Outputs are heatmaps saved out into pwd.
%
%


%% Input checks

if nargin ~= 4
    error('Function modularityRealVersusNull requires input args "realResFile", "nullResFile", "xlim" and "ylim"!');
end
if ~exist(realResFile, 'file') || ~exist(nullResFile, 'file')
    error('Input args "realResFile" and "nullResFile" must be valid file paths!');
end
if ~isvector(xlim) || numel(xlim)~=2
    error('Input arg "xlim" must be a numeric vector with two elements!');
end
if ~isvector(ylim) || numel(ylim)~=2
    error('Input arg "ylim" must be a numeric vector with two elements!');
end


%% Load data

realRes=load(realResFile);
nullRes=load(nullResFile);

% sanity checks - same methods and params?
if ~isequal(realRes.gammaValues, nullRes.gammaValues) || ~isequal(realRes.omegaValues, nullRes.omegaValues)
    error('Gamma or omega values are not the same across the two data sets!');
end
if ~strcmp(realRes.randmove, nullRes.randmove) || ~strcmp(realRes.method, nullRes.method)
    error('"Randmove" or "method" settings were not the same in the two data sets!');
end
if ~isequal(realRes.postprocess, nullRes.postprocess)
    error('"Postprocess" settings were not the same in the two data sets!');
end


%% Differences between real and null results

% get difference for main results
diffRes = struct;
fieldn = fieldnames(realRes.res);
for i = length(fieldn):-1:1
    if length(size(realRes.res.(fieldn{i}))) == 2  
        diffRes.(fieldn{i}) = realRes.res.(fieldn{i})-nullRes.res.(fieldn{i});
    else
        fieldn(i) = [];
    end
end


%% Heatmaps for differences

% axes ticks for heatmaps
xdata = realRes.omegaValues;
ydata = realRes.gammaValues;

% titles
titles = {'Mean no. of communities',...
    'SD of no. of communities',...
    'Mean quality (Q)',...
    'SD of quality (Q)',...    
    'Mean similarity (zRand score) across repetitions',...
    'SD of similarity (zRand score) across repetitions',...
    'Mean persistence',...
    'SD of persistence',...    
    };

% heatmaps for the differences
for i = 1:length(fieldn)
    
    f = figure;
    f.Position = [0 0 1800, 1080];
    heatmap(xdata(xlim(1):xlim(2)), ydata(ylim(1):ylim(2)), diffRes.(fieldn{i})(ylim(1):ylim(2), xlim(1):xlim(2)), 'Title', [fieldn{i}, ' - ', titles{i}], 'XLabel', 'omega', 'Ylabel', 'gamma');
    
    saveF = ['heatmap_modDifference_', fieldn{i},'_', num2str(xlim(1)), '_', num2str(xlim(2)), '_', num2str(ylim(1)), '_', num2str(ylim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end


%% Heatmaps for real results

for i = 1:length(fieldn)
    
    f = figure;
    f.Position = [0 0 1800, 1080];
    heatmap(xdata(xlim(1):xlim(2)), ydata(ylim(1):ylim(2)), realRes.res.(fieldn{i})(ylim(1):ylim(2), xlim(1):xlim(2)), 'Title', [fieldn{i}, ' - ', titles{i}], 'XLabel', 'omega', 'Ylabel', 'gamma');
    
    saveF = ['heatmap_modReal_', fieldn{i},'_', num2str(xlim(1)), '_', num2str(xlim(2)), '_', num2str(ylim(1)), '_', num2str(ylim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end


%% Heatmaps for null results

for i = 1:length(fieldn)
    
    f = figure;
    f.Position = [0 0 1800, 1080];
    heatmap(xdata(xlim(1):xlim(2)), ydata(ylim(1):ylim(2)), nullRes.res.(fieldn{i})(ylim(1):ylim(2), xlim(1):xlim(2)), 'Title', [fieldn{i}, ' - ', titles{i}], 'XLabel', 'omega', 'Ylabel', 'gamma');
    
    saveF = ['heatmap_modNull_', fieldn{i},'_', num2str(xlim(1)), '_', num2str(xlim(2)), '_', num2str(ylim(1)), '_', num2str(ylim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end



return

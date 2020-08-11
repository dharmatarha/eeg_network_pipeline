function diffRes = modularityRealVersusNull(realResFile, nullResFile, omegaLim, gammaLim)
%% Generating heatmaps from multiCommDetect(Null)Wrapper outputs
%
% USAGE: diffRes = modularityRealVersusNull(realResFile, 
%                                           nullResFile,
%                                           omegaLim=[1, length(omegaValues)], 
%                                           gammaLim=[1, length(gammmaValues)])
%
% Creates and saves a series of heatmaps to depict the results of
% corresponding multiCommDetectWrapper and multiCommDetectNullWrapper runs.
%
% 
% Mandatory inputs:
% realResFile   - Char array, path to file storing multiCommDetectWrapper
%               outputs
% nullResFile   - Char array, path to file storing multiCommDetectNullWrapper
%               outputs
%
% Optional inputs:
% omegaLim      - Numeric vector with two elements. X axis (omega axis) 
%               limits for the heatmaps. Defaults to full range of values 
%               (first and last elements of the omegaValues arrays from 
%               the multiCommDetectWrapper/multiCommDetectNullWrapper outputs).
% gammaLim      - Numeric vector with two elements. Y axis (gamma axis) 
%               limits for the heatmaps. Defaults to full range of values 
%               (first and last elements of the omegaValues arrays from 
%               the multiCommDetectWrapper/multiCommDetectNullWrapper outputs).
%
% Outputs are heatmaps saved out into pwd. There are three group of
% heatmaps:
% (1) Heatmaps for multiCommDetectWrapper outputs; saved into files
% starting with "heatmap_modReal_*"
% (2) Heatmaps for multiCommDetectNullWrapper outputs; saved into files
% starting with "heatmap_modNull_*"
% (3) Heatmaps for the difference between the two types of outputs; saved 
% into files starting with "heatmap_modDifference_*"
%
% The rest of the filenames describe the stat depicted (e.g. mean quality
% value) and the omega / gamma limits used (x /y axis limits). E.g.
% heatmap_modDifference_qMean_1_50_1_50.png
%


%% Input checks

if ~ismembertol(nargin, 2:4)
    error(['Function modularityRealVersusNull requires input args "realResFile" ',...
        'and "nullResFile" while args "omegaLim" and "gammaLim" are optional!']);
end
if ~exist(realResFile, 'file') || ~exist(nullResFile, 'file')
    error('Input args "realResFile" and "nullResFile" must be valid file paths!');
end
if nargin == 2
    omegaLim = false;
    gammaLim = false;
elseif nargin == 3
    omegaLim = false;
elseif nargin == 4
    if ~isvector(omegaLim) || numel(omegaLim)~=2
        error('Input arg "omegaLim" must be a numeric vector with two elements!');
    end
    if ~isvector(gammaLim) || numel(gammaLim)~=2
        error('Input arg "gammaLim" must be a numeric vector with two elements!');
    end
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

% get omega and gamma limits if they were not supplied
if ~omegaLim
    omegaLim = [1, length(realRes.omegaValues)];
end
if ~gammaLim
    gammaLim = [1, length(realRes.gammaValues)];
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
    heatmap(xdata(omegaLim(1):omegaLim(2)),... 
        ydata(gammaLim(1):gammaLim(2)),... 
        diffRes.(fieldn{i})(gammaLim(1):gammaLim(2), omegaLim(1):omegaLim(2)),... 
        'Title', [fieldn{i}, ' - ', titles{i}],... 
        'XLabel', 'omega',... 
        'Ylabel', 'gamma',...
        'ColorMap', jet);
    
    saveF = ['heatmap_modDifference_', fieldn{i},'_', num2str(omegaLim(1)), '_', num2str(omegaLim(2)), '_', num2str(gammaLim(1)), '_', num2str(gammaLim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end


%% Heatmaps for real results

for i = 1:length(fieldn)
    
    f = figure;
    f.Position = [0 0 2400, 1400];
    heatmap(xdata(omegaLim(1):omegaLim(2)),... 
        ydata(gammaLim(1):gammaLim(2)),... 
        realRes.res.(fieldn{i})(gammaLim(1):gammaLim(2), omegaLim(1):omegaLim(2)),... 
        'Title', [fieldn{i}, ' - ', titles{i}],... 
        'XLabel', 'omega',... 
        'Ylabel', 'gamma',...
        'ColorMap', jet);
    
    saveF = ['heatmap_modReal_', fieldn{i},'_', num2str(omegaLim(1)), '_', num2str(omegaLim(2)), '_', num2str(gammaLim(1)), '_', num2str(gammaLim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end


%% Heatmaps for null results

for i = 1:length(fieldn)
    
    f = figure;
    f.Position = [0 0 1800, 1080];
    heatmap(xdata(omegaLim(1):omegaLim(2)),... 
        ydata(gammaLim(1):gammaLim(2)),... 
        nullRes.res.(fieldn{i})(gammaLim(1):gammaLim(2), omegaLim(1):omegaLim(2)),... 
        'Title', [fieldn{i}, ' - ', titles{i}],... 
        'XLabel', 'omega',... 
        'Ylabel', 'gamma',...
        'ColorMap', jet);
    
    saveF = ['heatmap_modNull_', fieldn{i},'_', num2str(omegaLim(1)), '_', num2str(omegaLim(2)), '_', num2str(gammaLim(1)), '_', num2str(gammaLim(2)),'.png'];
    
    saveas(f, saveF);
    
    close(f);
    
end



return

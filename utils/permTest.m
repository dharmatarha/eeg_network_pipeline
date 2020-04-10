function [pEst, realDiff, permDiff, CohenD] = permTest(a, b, varargin)
%% Permutation test for the difference between two data vectors
%
% USAGE: [pEst, realDiff, permDiff, CohenD] = permTest(a, b, perm = 10000, stat = 'mean')
%
% Simple permutation test to estimate the statistical significance of the
% difference between two data sets in terms of the mean or the median.
% Provides the one-sided estimate, no correction for the double test.
% The script relies on permuting data indices.
%
% Optional input args "perm" and "stat" are inferred from the data types of
% the inputs.
% 
% Mandatory inputs:
% a, b      - Data vectors.
%
% Optional inputs:
% perm      - No. of permutations used for estimation. Defaults to 10^4.
% stat      - Test statistic to use, string. One of {'mean', 'median',
%           'std'}
%
% Outputs:
% pEst      - Estimated probability of the null hypothesis (i.e. that
%           stat(a) = stat(b)).
% realDiff  - Real difference between data "a" and "b" in terms of the test
%           statistic.
% permDiff  - Differences from random permutations.
% CohenD    - Cohen's d (effect size). Note that Cohen's d relies on the
%           mean irrespective of the stat used for the permutation test
%
%

%% Input checks

% need at least two args
if nargin < 2
    error('Function permTest requires at least two input args, data vectors "a" and "b"!');
end
% check optional args
if ~isempty(varargin)
    % infer which one the first element is
    if ischar(varargin{1})
        stat = varargin{1};
    elseif isnumeric(varargin{1})
        perm = varargin{1};
    end
    % check if ther is another optional input
    if length(varargin) == 2
        if exist('stat', 'var') && isnumeric(varargin{2})
            perm = varargin{2};
        elseif exist('perm', 'var') && ischar(varargin{2})
            stat = varargin{2};
        else 
            error(['Could not parse second optional argument!',... 
                'Input arg "stat" should be a string, input arg "perm" should be numeric!']);
        end
    elseif length(varargin) > 2
        error('Too many optional input args, cannot parse! Optional args are "stat" and "perm"!');
    end
end
% check optionals
if ~exist('stat', 'var')
    stat = 'mean';
else
    if ~ismember(stat, {'mean', 'median', 'std'})
        error('Input arg "stat" should be one of "mean", "median" or "std"!');
    end
end
if ~exist('perm', 'var')
    perm = 10000;
else
    if ~ismembertol(perm, 1:1000000)
        error('Input arg "perm" should be between 1 and 10^6!');
    end    
end
    
% chekc mandatories    
if ~isvector(a) || ~isvector(b)
    error('Input args "a" and "b" should be numeric vectors!');
end
if ~iscolumn(a)
    a = a';
end
if ~iscolumn(b)
    b = b';
end

% user message
disp([char(10), 'Called permTest for data vectors sized: ',...
    char(10), num2str(size(a)), ...
    char(10), num2str(size(b)),...
    char(10), 'No. of permutations: ', num2str(perm),...
    char(10), 'Test statistic: ', stat]);


%% permute data, calculate test stat

% check for NaN
if any(isnan(a)) || any(isnan(b))
    warning(['There is at least one NaN value in the data vectors. ',...
        'We treat this with ''omitnan'' flags but you should know about this!']);
end

% real difference
switch stat
    case 'mean'
        realDiff = mean(a, 'omitnan')-mean(b, 'omitnan');
    case 'median'
        realDiff = median(a, 'omitnan')-median(b, 'omitnan');
    case 'std'
        realDiff = std(a, 'omitnan')-std(b, 'omitnan');   
end

% prepare for permutations
data = [a; b];  % concatenate data vectors
la = length(a);  % sample no for vector "a"
lb = length(b);  % sample no for vector "b"
ld = la + lb;  % no. of all samples
% preallocate results var
permDiff = nan(perm, 1);

% permutations loop
for i = 1:perm
    % permute
    d = data(randperm(ld));
    % get test stat for permuted data
    switch stat
        case 'mean'
            permDiff(i) = mean(d(1:la), 'omitnan')-mean(d(la+1:ld), 'omitnan');
        case 'median'
            permDiff(i) = median(d(1:la), 'omitnan')-median(d(la+1:ld), 'omitnan');   
        case 'std'
            permDiff(i) = std(d(1:la), 'omitnan')-std(d(la+1:ld), 'omitnan');
    end
end


%% get p value and effect size

if realDiff <= 0
    pEst = 1-(sum(permDiff>realDiff)/perm);
elseif realDiff > 0
    pEst = 1-(sum(permDiff<realDiff)/perm);
end

CohenD = (mean(a)-mean(b))/(((std(a)^2+std(b)^2)/2)^0.5);

% user message
disp([char(10), 'Real difference was ', num2str(realDiff),...
    char(10), 'Estimated probability of H0 (equality): ', num2str(pEst),...
    char(10), 'Effect size (Cohen''s d): ', num2str(CohenD), char(10)]);


return



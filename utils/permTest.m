function [pEst, realDiff, permDiff, CohenD, studentDiff] = permTest(a, b, varargin)
%% Permutation test for the difference between two data vectors
%
% USAGE: [pEst, realDiff, permDiff, CohenD] = permTest(a, b, 
%                                                       perm = 10000, 
%                                                       stat='mean', 
%                                                       student='studentized',
%                                                       verbosity='verbose')
%
% Two sample random permutation test to estimate the statistical 
% significance of the difference between two data sets in terms of the 
% mean, SD or the median.
% Provides the one-sided estimate, no correction for the double test.
% The script relies on permuting data indices.
%
% Optional input args "perm" and "stat" are inferred from the data types of
% the inputs.
% 
% Test statistic (the sample function outcome under consideration) is
% studentized by default (see the literature at the end of the help for
% more details). 
% IMPORTANT: STUDENTIZATION ONLY WORKS NOW FOR THE
% DIFFERENCE OF MEANS AS TEST STATISTIC.
%
% Mandatory inputs:
% a, b      - Data vectors.
%
% Optional inputs:
% perm      - Numeric value, no. of permutations used for estimation. In 
%           range 1:10^6. Defaults to 10^4.
% stat      - Char array, one of {'mean', 'median', 'std'}. Test 
%           statistic to use. Defaults to 'mean'.
% student   - Char array, one of {'studentized', 'raw'}. Controls if test
%           statistic is studentized or not. Defaults to 'studentized'.
% verbosity - Char array, one of {'verbose', 'silent'}. Controls verbosity,
%           i.e., whether user messages are displayed or the script runs 
%           silently. Defaults to 'verbose'.
%
% Outputs:
% pEst      - Estimated probability of the null hypothesis (i.e. that
%           stat(a) = stat(b)).
% realDiff  - Real difference between data "a" and "b" in terms of the test
%           statistic.
% permDiff  - Differences from random permutations.
% CohenD    - Cohen's d (effect size). Note that Cohen's d relies on the
%           mean irrespective of the stat used for the permutation test
% studentDiff - Studentized test statistic for the real difference. Empty
%           if "student" was "raw".
%
%
% NOTES:
% Random permutation tests are not liked universally. For example, Rand
% Wilcox in his "Modern Statistics for the Social and Behavioral Sciences.
% A Practical Introduction (2nd ed.)" argues against using them, because 
% they are not robust/powerful enough when comparing different /
% heavy-tailed distributions. However, a simple studentization does wonders
% (asymptotically), see e.g.:
% (1) Chung, E., & Romano, J. P. (2013). Exact and asymptotically 
%       robust permutation tests. The Annals of Statistics, 41(2), 484-507.
% (3) Janssen, A. (1997). Studentized permutation tests for non-iid 
%       hypotheses and the generalized Behrens-Fisher problem. 
%       Statistics & probability letters, 36(1), 9-21.
%

%% Input checks

% need at least two args
if ~ismember(nargin, 2:6)
    error(['Function permTest requires input args "a" and "b" (data vectors) ',...
        'while input args "perm", "stat", "student" and "verbosity" are optional (see the help)!']);
end
% check mandatory args    
if ~isvector(a) || ~isvector(b)
    error('Input args "a" and "b" should be numeric vectors!');
end
if ~iscolumn(a)
    a = a';
end
if ~iscolumn(b)
    b = b';
end
% check optional args
if ~isempty(varargin)
    % loop through each varargin element
    for v = 1:length(varargin)
        % sort by type and range
        if ischar(varargin{v}) && ismember(varargin{v}, {'mean', 'median', 'std'}) && ~exist('stat', 'var')
            stat = varargin{v};
        elseif isnumeric(varargin{v}) && ismembertol(varargin{v}, 1:1000000) && ~exist('perm', 'var')
            perm = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'studentized', 'raw'}) && ~exist('student', 'var')
            student = varargin{v};              
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'verbose', 'silent'}) && ~exist('verbosity', 'var')
            verbosity = varargin{v};    
        else 
            error('At least one optional input arg could not be mapped nicely to "perm", "stat", "student" or "verbosity"!');
        end
    end
end
% assign default values where necessary
if ~exist('stat', 'var')
    stat = 'mean';
end
if ~exist('perm', 'var')
    perm = 10000;  
end
if ~exist('student', 'var')
    student = 'studentized';  
end
if ~exist('verbosity', 'var')
    verbosity = 'verbose';  
end
% turn verbosity into logical
if strcmp(verbosity, 'verbose')
    verbosity = true;
else
    verbosity = false;
end

% user message
if verbosity
    disp([char(10), 'Called permTest for data vectors sized: ',...
        char(10), num2str(size(a)), ...
        char(10), num2str(size(b)),...
        char(10), 'No. of permutations: ', num2str(perm),...
        char(10), 'Test statistic: ', stat,...
        char(10), 'Studentized test statistic: ', student]);
end

% turn student into logical
if strcmp(student, 'studentized')
    student = true;
else
    student = false;
end


%% permute data, calculate test stat

% check for NaN
if verbosity
    if any(isnan(a)) || any(isnan(b))
        warning(['There is at least one NaN value in the data vectors. ',...
            'We treat this with ''omitnan'' flags but you should know about this!']);
    end
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

% studentized versions of differences
% ONLY MEAN IS SUPPORTED THIS WAY AS OF NOW
if student
    switch stat
        case 'mean'
            studentDiff = ld^0.5*(mean(a, 'omitnan')-mean(b, 'omitnan'))/...
                (sqrt(ld*(var(a)/la+var(b)/lb)));
    end
else
    studentDiff = [];
end

% permutations loop
for i = 1:perm
    
    % permute
    d = data(randperm(ld));
    
    % get test stat for permuted data
    
    % if studentized test statistic is to be used
    if student
        pa = d(1:la);
        pb = d(la+1:ld);        
        switch stat
            case 'mean'
                permDiff(i) = ld^0.5*(mean(pa, 'omitnan')-mean(pb, 'omitnan'))/...
                    (sqrt(ld*(var(pa)/la+var(pb)/lb)));
        end 
        
    % if non-studentized, "raw" statistic is used
    else
        switch stat
            case 'mean'
                permDiff(i) = mean(d(1:la), 'omitnan')-mean(d(la+1:ld), 'omitnan');
            case 'median'
                permDiff(i) = median(d(1:la), 'omitnan')-median(d(la+1:ld), 'omitnan');   
            case 'std'
                permDiff(i) = std(d(1:la), 'omitnan')-std(d(la+1:ld), 'omitnan');
        end  % switch stat
        
    end  % if student
    
end  % for i


%% get p value and effect size

if student
    if studentDiff <= 0
        pEst = 1-(sum(permDiff>studentDiff)/perm);
    elseif realDiff > 0
        pEst = 1-(sum(permDiff<studentDiff)/perm);
    end    
else    
    if realDiff <= 0
        pEst = 1-(sum(permDiff>realDiff)/perm);
    elseif realDiff > 0
        pEst = 1-(sum(permDiff<realDiff)/perm);
    end
end

CohenD = (mean(a)-mean(b))/(((std(a)^2+std(b)^2)/2)^0.5);

% user message
if verbosity
    disp([char(10), 'Real difference was ', num2str(realDiff),...
        char(10), 'Estimated probability of H0 (equality): ', num2str(pEst),...
        char(10), 'Effect size (Cohen''s d): ', num2str(CohenD), char(10)]);
end


return



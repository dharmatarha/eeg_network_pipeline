function [h, pCrit] = fdr(pValues, varargin)
%% False-Discovery Rate calculation
%
% USAGE: [h, pCrit] = fdr(pValues, q = 0.05, method = 'bh')
%
% Calculates false discovery rate for supplied p (significance) values and
% supplied q. Method might be Benjamini–Hochberg ('bh') or 
% Benjamini–Yekutieli ('by').
%
% Mandatory input:
% pValues       - Vector of uncorrected significance values.
%
% Optional inputs:
% q             - False discovery rate to control for. Defaults to 0.05
% method        - Procedure for FDR. One of {'bh', 'by'}, referring to
%               Benjamini-Hochberg and Benjamini-Yekutieli, respectively.
%               Defaults to 'bh'. Benjamini-Hochberg is valid under the
%               independence or positive dependence assumption,
%               Benjamini-Yekutieli is always valid but less powerful.
%
% Outputs:
% h             - Null hypothesis rejection values for all element in
%               pValues. If 1, the null hypothesis for that test is
%               rejected, i.e. the finding is significant.
% pCrit         - Critical p value. Any pValues<pCrit is significant.
%
%


%% Input checks

% check for mandatory input
if nargin < 1 
    error('Function fdr needs input arg "pValues" and optional input args "q" and "method"!');
end

% sort optional args into "q" and "method"
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ~exist('method', 'var')
            method = varargin{v};
        elseif isnumeric(varargin{v}) && ~exist('q', 'var')
            q = varargin{v};
        else
            error('Could not parse at least one optional input arg!');
        end
    end
end

% assign default values
if ~exist('method', 'var')
    method = 'bh';
end
if ~exist('q', 'var')
    q = 0.05;
end

% check values
if q<=0 || q>=1
    error(['Input arg "q" = ', num2str(q), ', looks very fishy!']);
end
if ~ismember(method, {'bh', 'by'})
    error(['Input arg "method" should be one of {''bh'', ''by''}, not ', method, '!']);
end
if ~isvector(pValues)
    error('Input arg "pValues" should be a numerical vector!');
end
if iscolumn(pValues)
    pValues = pValues';
end


%% FDR procedure

% sort p values
pSorted = sort(pValues);
% no. of tests
m = length(pSorted);

% Benjamini-Hochberg, for independent or positively dependent tests
if strcmp(method, 'bh')
    % get threshold values for all sorted ps
    thValues=(1:m)*q/m;
    
% Benjamini-Yekutieli, for arbitrary dependence structure
elseif strcmp(method, 'by')
    % there is an extra correction term in the denominator: harmonic series up to m
    harmonicM = sum(1./(1:m));
    thValues=(1:m)*q/(m*harmonicM);
    
end

% compare the threshold values to the significance values, look for
% largest passing the test
pSortedLargest = find(pSorted <= thValues, 1, 'last');

if ~isempty(pSortedLargest)
    % critical p
    pCrit = pSorted(pSortedLargest);
    % hypothesis tests
    h = (pValues <= pCrit);
else
    % otherwise return zeros 
    pCrit = 0;
    h = zeros(1, length(pValues));
end


return


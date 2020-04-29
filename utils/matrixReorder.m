function [mNew, old2new] = matrixReorder(m, oldLabels, newLabels)
%% Reordering a square matrix based on a new ordering of row/column labels
%
% USAGE [mNew, old2new] = matrixReorder(m, oldLabels, newLabels)
%
% Helper function for reordering connectivity / adjacency matrices on the
% basis of ROI lists (e.g. anatomical regions). 
%
% Inputs:
% m         - Numerical square matrix.
% oldLabels - Cell array of row (and column) labels for input matrix "m"
% newLabels - Cell array of row (and column) labels of output matrix "mNew"
%
% Outputs:
% mNew      - Numerical square matrix, reordered version of input arg "m".
%           I.e., m(old2new, old2new) == mNew
% old2new   - Numerical vector, indices for transforming input arg 
%           "oldLabels" to "newLabels", so that 
%           oldLabels(old2new) == newLabels 
%


%% Input checks

if nargin ~= 3
    error('Function matrixReorder needs input args "m", "oldLabels" and "newLabels"!');
end
if ~ismatrix(m) || ~isequal(size(m, 1), size(m, 2))
    error('Input arg "m" should be a square matrix!');
end
if (~isvector(oldLabels) && ~iscell(oldLabels)) || (~isvector(newLabels) && ~iscell(newLabels))
    error('Input args "oldLabels" and "newLabels" should be cell vectors!');
end
if ~isequal(size(oldLabels), size(newLabels))
    error('Input args "oldLabels" and "newLabels" must have same size!');
end


%% Compare label cell arrays, get indices from the old to the new

% make sure label args are both column vectors
if ~iscolumn(oldLabels)
    oldLabels = oldLabels';
end
if ~iscolumn(newLabels)
    newLabels = newLabels';
end

% get indices transforming old to new
[c, ~, old2new] = intersect(newLabels, oldLabels, 'stable');

% sanity check - are the elements in the two label arrays the same?
if ~isequal(length(c), length(oldLabels))
    error('Input args "oldLabels" and "newLabels" have at least one different element!');
end


%% Re-order matrix

mNew = m(old2new, old2new);


return


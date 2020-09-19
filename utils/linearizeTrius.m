function dataLin = linearizeTrius(dataArray, k)
%% Helper function to vectorize (linearize) the upper triangles of a set of matrices.
%
% USAGE: dataLin = linearizeTrius(dataArray, k=1)
%
% The standard use case is to linearize connectivity matrices (across 
% epochs and conditions) so that we can use  simple vector-based measures /
% transformations on them (e.g. correlation).
%
% Mandatory input:
% dataArray         - 3D or 4D numeric array, with the first two dimensions
%               having equal sizes (square matrices).  
%
% Optional input:
% k                 - Numeric value, the "K" input arg supplied to the
%               built-in function "triu". Determines the diagonal for the
%               upper triangle. 
%               From "help triu": triu(X,K) is the 
%               elements on and above the K-th diagonal of X.  
%               K = 0 is the main diagonal, K > 0 is above the main
%               diagonal and K < 0 is below the main diagonal. 
%               Defaults to 1.
% 
% Output:
% dataLin           - 2D or 3D numeric array. The function linearizes the 
%               first two dimensions of the input array "dataArray", so the
%               output "dataLin" has one less dimension than "dataArray".
%


%% Input checks

if ~ismember(nargin, 1:2)
    error('Function linTrius requires input arg "dataArray", while input arg "k" is optional!');
end
if ~isnumeric(dataArray) || ~ismember(length(size(dataArray)), 3:4)
    error('Input arg "dataArray" should be a 3D or 4D numeric array!');
end
if nargin == 1
    k = 1;
else
    if ~isnumeric(k) || mod(k, 1)~=0 || numel(k)~=1
        error('Input arg "k" must be an integer!');
    end
end
    

%% Get size of input data + preallocate output var

% size of input data
if length(size(dataArray)) == 3
    dimNo = 3;
    [l1, l2, l3] = size(dataArray);
else 
    dimNo = 4;
    [l1, l2, l3, l4] = size(dataArray);
end

% check if "k" is sensible given the input array size
if k >= l2
    error(['Input arg "k" is equal to or larger than the second dim of the ',...
        'input array, upper triangle of "k"-th diagonal will be empty!']);
end

% no. of elements in upper triangles of matrices defined by l1 and l2, with
% k-th diagonal
linNo = sum(sum(triu(ones(l1, l2), k)));
% preallocate output var
if dimNo==3
    dataLin = zeros(linNo, l3);
else
    dataLin = zeros(linNo, l3, l4);
end


%% Extract upper triangles in for loops 

if dimNo == 3
    for i = 1:l3
        % get upper triangular part
        tmp = dataArray(:, :, i);
        idx = triu(true(size(tmp)), k); 
        dataLin(:, i) = tmp(idx);
    end  % for i
else
    for j = 1:l4
        for i = 1:l3
            % get upper triangular part
            tmp = dataArray(:, :, i, j);
            idx = triu(true(size(tmp)), k); 
            dataLin(:, i, j) = tmp(idx);
        end  % for i
    end  % for j
end  % if dimFlag


return
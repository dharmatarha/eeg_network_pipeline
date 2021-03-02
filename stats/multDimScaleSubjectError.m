function [subjectErrors] = multDimScaleSubjectError(distanceMatrix, subjectCoordinateMatrix, varargin)
%% Testing the error of multidimensional scaling for each subject
%
% USAGE: subjectErrors = multDimScaleSubjectError(distanceMatrix, subjectCoordinateMatrix, numberOfDimensions=size(subjectCoordinateMatrix, 2), verbose = true)
%
% If given a matrix of distances between subjects and the matrix of
% estimated coordinates for each subject (result of multidimensional 
% scaling), the function calculates the error of reconstructed Euclidean
% distances for each subject.
%
% Workflow:
% (1) Calculate reconstructed Euclidean distances between subjects
% based on the coordinates estimated by multidimensional scaling.
% (2) Calculate the difference between the original distance matrix
% and the matrix of reconstructed distances
% (3) Sum the absolute differences (errors) for each subject
%
% Mandatory inputs:
% distanceMatrix           - Numeric square matrix with dimensions: subjects X subjects.
%                          Contains distances between subjects. Must be symmetric, values
%                          on the diagonal must be zeros.
% subjectCoordinateMatrix  - Numeric matrix, matrix of coordinates estimated by 
%                          multidimensional scaling. The number of rows must be equal to
%                          that of distanceMatrix.
%
% Optional inputs:
% numberOfDimensions       - Number of dimensions to consider for distance reconstruction.
%                          Defaults to all dimensions.
% verbose                  - Logical value (true or false). Verbosity, "false" 
%                          meaning no user messages, "true" meaning user messages.
%
% Output:
% subjectErrors            - Row vector, one numeric value for each subject. Contains
%                          the summed error of reconstructed Euclidean distances for each subject,
%                          devided by the original distances.
%


%% Input checks

% check number of arguments
if ~ismember(nargin, 2:4)
    error(['Function multDimScaleSubjectError requires input args "distanceMatrix" and "subjectCoordinateMatrix", ',...
        'while input args "numberOfDimensions" and "verbose" are optional!']);
end
% check mandatory args
if ~isnumeric(distanceMatrix) || ~ismatrix(distanceMatrix) || size(distanceMatrix, 1)~=size(distanceMatrix, 2)
    error(['Input arg "distanceMatrix" should be a numeric square matrix ',...
        '(matrix of distances between subjects)!']);
end
if ~isnumeric(subjectCoordinateMatrix) || ~ismatrix(subjectCoordinateMatrix) || size(subjectCoordinateMatrix, 1)~=size(distanceMatrix, 1)
    error(['Input arg "subjectCoordinateMatrix" should be a numeric matrix ',...
        'with the same number of rows as "distanceMatrix" (matrix of coordinates estimated by multidimensional scaling)!']);
end

% check optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ~exist('numberOfDimensions', 'var')
            numberOfDimensions = varargin{v};
        elseif islogical(varargin{v}) && numel(varargin{v})==1 && ~exist('verbose', 'var')
            verbose = varargin{v};
        else
            error('An optional input arg does not match nicely to "numberOfDimensions" or "verbose"!');
        end
    end
end
% assign defaults
if ~exist('numberOfDimensions', 'var')
    numberOfDimensions = size(subjectCoordinateMatrix, 2);
end
if ~exist('verbose', 'var')
    verbose = true;
end
% any further check
if any(diag(distanceMatrix))
    error('There is at least one non-zero value on the diagonal of "distanceMatrix"!');
end
if numberOfDimensions > size(subjectCoordinateMatrix, 2)
    numberOfDimensions = size(subjectCoordinateMatrix, 2);
    warning('numberOfDimensions limited to the second dimension of the subjectCoordinateMatrix');
end

% user message if verbose
if verbose
    disp([char(10), 'Called multDimScaleSubjectError function with input args: ',...
        char(10), 'Distance matrix of size ', num2str(size(distanceMatrix)),...
        char(10), 'numberOfDimensions (number of dimensions to consider): ', num2str(numberOfDimensions),...
        char(10), 'Verbosity: ', num2str(verbose), char(10)]);
end


%% Calculate distance error for each subject

reconstructedDistanceMatrix = squareform(pdist(subjectCoordinateMatrix(:, 1:numberOfDimensions)));

distanceErrorMatrix = distanceMatrix - reconstructedDistanceMatrix;

% normalizedDistanceErrorMatrix = abs(distanceErrorMatrix ./ distanceMatrix);
% normalizedDistanceErrorMatrix(isnan(normalizedDistanceErrorMatrix)) = 0;
% subjectErrors = sum(normalizedDistanceErrorMatrix);

subjectErrors = sum(abs(distanceErrorMatrix)) ./ sum(abs(distanceMatrix));

end
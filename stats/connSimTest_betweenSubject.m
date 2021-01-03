function simRes = connSimTest_betweenSubject(connArray, varargin)
%% Testing the similarity of connectivity data between different subjects
%
% USAGE: simRes = connSimTest_betweenSubject(connArray, metric='corr') 
%
% If given connectivity (adjacency) matrices for a set of epochs, across
% multiple subjects (in input arg "connArray"), the function calculates the
% consistency of average connectivity matrices between different subjects.
%
% Workflow:
% (1) Calculate average connectivity (adjacency) matrices for each subject
% (averaging is done across all epochs of a given subject)
% (2) Calculate the similarity of two averaged conenctivity matrices
% for all possible pairing of all subjects.
%
% At the moment, supports similarity / distance metrics (1) correlation 
% (input arg "metric" = 'corr'), (2) Eucledian distance ('eucl'), and (3)
% DeltaCon ('deltaCon', see /networkSimilarity/deltaCon.m for details).
%
% Mandatory inputs:
% connArray     - 4D numeric array, with dimensions: subjects X epochs X
%               ROIs X ROIs. Contains a connectivity matrix for each epoch
%               of each subject. Might only contain valid data in the upper
%               triangles of the connectivity matrices, the rest might be
%               NaN (in case of symmetric connectivity measures).
%
% Optional inputs:
% metric        - Char array, one of {'corr', 'eucl', 'deltaCon'}.
%               Similarity metric for comparing connectivity matrices.
%               DeltaCon relies on the similarly named function in
%               /networkSimilarity. Defaults to 'corr'.
%
% Output:
% simRes        - 2D numeric array sized subjects X subjects. Contains
%               connectivity matrix similarity values for all possible
%               pairings of all subjects. Only upper triangle is populated.
%


%% Input checks

% cehck no. if inputs
if ~ismember(nargin, 1:2)
    error('Function connSimTest_betweenSubject requires input arg "connArray" while the arg "metric" is optional!');
end
% check mandatory input
if ~isnumeric(connArray) || numel(size(connArray))~=4 || size(connArray, 3)~=size(connArray ,4)
    error('Input arg "connArray" should be a 4D numeric array where the 3rd and 4th dimensions have the same size!');
end
% check optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v}) && ismember(varargin{v}, {'corr', 'eucl', 'deltaCon'}) && ~exist('metric', 'var')
            metric = varargin{v};
        else
            error('Input could not be mapped nicely to arg "metric"!');
        end
    end
end
% assign default values if necessary
if ~exist('metric', 'var')
    metric = 'corr';
end

% user message
disp([char(10), 'Called connSimTest_betweenSubject with input args: ',...
    char(10), 'Input array sized ', num2str(size(connArray)), ...
    char(10), 'Similarity metric: ', metric]);


%% Calculate similarities

% clock for whole run
funcClock = tic;

% get subject, epoch and channel/ROI numbers
[subNo, ~, roiNo, ~] = size(connArray);

% preallocate results var
simRes = nan(subNo, subNo);

% average over epochs before subject-to-subject
% comparisons
connArray = squeeze(mean(connArray, 2));

% user message
disp([char(10), 'Calculating...']);

% loop through subjects
for subIdx = 1:subNo
    
    % subject-level clock
    subClock = tic;
    
    % select subject's connectivity data (averaged over epochs)
    subMatrixA = squeeze(connArray(subIdx, :, :));
    
    % loop through compared subjects
    for compSubIdx = subIdx+1 : subNo
        
        % select compared subject's connectivity data (averaged over
        % epochs)
        subMatrixB = squeeze(connArray(compSubIdx, :, :));
        
        % for certain metrics get symmetric adjacency matrices with zeros at diagonal
        if ismember(metric, {'eucl', 'deltaCon'})
            subMatrixA = triu(subMatrixA, 1) + triu(subMatrixA, 1)'; 
            subMatrixB = triu(subMatrixB, 1) + triu(subMatrixB, 1)';
            % also standardize the scale of connections across the two
            % matrices to a common sum (=10)
            subMatrixA = 10*subMatrixA./sum(subMatrixA(:));
            subMatrixB = 10*subMatrixB./sum(subMatrixB(:));            
        end

        % calculate similarity according to arg "metric"
        switch metric
            
            case 'corr'
                % linearize upper triangles of mean connectivity matrices
                linA = subMatrixA(triu(true(roiNo), 1));
                linB = subMatrixB(triu(true(roiNo), 1));
                simRes(subIdx, compSubIdx) = corr(linA, linB);              
                
            case 'eucl'
%                 % similarity is based on norm of difference
%                 simRes(subIdx, compSubIdx) = 1/(1 + norm(subMatrixA-subMatrixB, 'fro'));
                % we use a distance measure if 'eucl' is selected
                % (Frobenius norm)
                simRes(subIdx, compSubIdx) = norm(subMatrixA-subMatrixB, 'fro');                  
                
            case 'deltaCon'
%                 % first output arg of deltaCon is similarity
%                 simRes(subIdx, compSubIdx) = deltaCon(subMatrixA, subMatrixB, false);  % "false" is for verbosity
                % we use a distance measure (DeltaCon distance) if
                % 'deltaCon' is selected (second output of deltaCon.m)
                [~, simRes(subIdx, compSubIdx)] = deltaCon(subMatrixA, subMatrixB, false);  % "false" is for verbosity                
                
        end  % switch metric
        
    end  % for permIdx
    
    % user message
    disp([char(10), 'Done with subject ', num2str(subIdx),... 
        ', took ', num2str(round(toc(subClock), 2)), ' secs']);
    
end  % for subIdx

%% End message, return

% user message
disp([char(10), 'Done with everything, running the function took ',... 
    num2str(round(toc(funcClock), 2)), ' secs']);


return
        
        













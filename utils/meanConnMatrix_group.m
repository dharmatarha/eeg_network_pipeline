function meanConnMatrix = meanConnMatrix_group(connArray)
%% Calculating the mean connectivity matrix on the group level
%
% USAGE: simRes = meanConnMatrix_group(connArray) 
%
% If given connectivity (adjacency) matrices for a set of epochs, across
% multiple subjects (in input arg "connArray"), the function calculates the
% average connectivity matrix on the group level.
%
% Input:
% connArray     - 4D numeric array, with dimensions: subjects X epochs X
%               ROIs X ROIs. Contains a connectivity matrix for each epoch
%               of each subject. Might only contain valid data in the upper
%               triangles of the connectivity matrices, the rest might be
%               NaN (in case of symmetric connectivity measures).
%
% Output:
% meanConnMatrix - 2D numeric array sized (number of ROIs) X (number of ROIs).
%


%% Input checks

% check mandatory input
if ~isnumeric(connArray) || numel(size(connArray))~=4 || size(connArray, 3)~=size(connArray ,4)
    error('Input rag "connArray" should be a 4D numeric array where the 3rd and 4th dimensions have the same size!');
end


%% Calculation

% clock for whole run
funcClock = tic;

% average over epochs before subject-to-subject
% comparisons
connArray = squeeze(mean(connArray, 2));

% user message
disp([char(10), 'Calculating...']);

% average subjects' data
meanConnMatrix = squeeze(mean(connArray, 1));


%% End message, return

% user message
disp([char(10), 'Done with everything, running the function took ',... 
    num2str(round(toc(funcClock), 2)), ' secs']);


return
        
        













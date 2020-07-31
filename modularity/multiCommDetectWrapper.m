function multiCommDetectWrapper(realConn, nullConn, gammaValues, omegaValues, varargin)
%% Wrapper for calling genlouvain on a multilayer connectivity matrix repeatedly
%
% USAGE: multiCommDetectWrapper(realConn, 
%                               nullConn, 
%                               gammaValues, 
%                               omegaValues, 
%                               rep=100, 
%                               method='iterated', 
%                               randmove='moverandw', 
%                               postprocess='postprocess-ordinal-multilayer')
%
% The function explores a range of resolution parameters (spatial and
% temporal, that is, gamma and omega, respectively) for a multilayer
% (usually temporal) connectivity data set. It calls the functions 
% getMultiLayerConnMatrix.m, genlouvain.m (from GenLouvain toolbox), 
% consensus_similarity.m and zrand.m (both from Network Community Toolbox) 
% repeatedly.
%
% The function does not save or return all the partitions generated.
% Instead it calculates a number of statistics for each gamma-omega
% pairing: mean and std of no. of communities; mean and std of Q; estimates
% for mean and std of partition-distance (z-rand score); mean and std of
% persistence;
% It also calculates and returns the consensus similarity partitions (one
% partition per (gamma, omega) pair)
% 
% Mandatory inputs:
% realConn    - 3D numeric array, contains a set of connectivity matrices.
%               Its dimensions are (nodeNo, nodeNo, epochNo).
% nullConn    - 2D or 3D numeric array, contains null (e.g. surrogate) 
%               connectivity matrix / matrices. If 2D, the same null 
%               matrix is used for all real connectivity matrices.
%               Its dimensions are (nodeNo, nodeNo (, epochNo)).
% gammaValues - Numeric vector, range of spatial resolution parameters to 
%               explore (e.g. 0.1:0.1:2).
% omegaValues - Numeric vector, range of temporal resolution parameters 
%               (inter-layer edge weights) to explore 
%               (e.g. logspace(-4, -3, 100)). Omega determines uniform
%               within-node inter-layer connections (all inter-layer
%               edges have the same weight). 
% rep         - Numeric value, number of repetitions for running
%               genlouvain.m with given gamma and omega. Must be a positive
%               integer in the range [10:10^4];
%


%% Input checks

% overall number of args
if nargin ~= 5
    error(['Function multiCommDetectWrapper requires input args "realConN", '... 
        '"nullConn", "gammaValues", omegaValues" and "rep"!']);
end
% check each arg
if ~isnumeric(realConn) || length(size(realConn)) ~= 3 || size(realConn, 1) ~= size(realConn, 2)
    error(['Input arg "realConN" should be a 3D numeric array with the first ',...
        'two dimensions having euqal size (nodeNo*nodeNo*epochNo)!']);
end
if ~isnumeric(nullConn) || ~ismember(length(size(nullConn)), 2:3)  || size(nullConn, 1) ~= size(nullConn, 2)
    error(['Input arg "nullConN" should be a 2D or 3D numeric array with the first ',...
        'two dimensions having euqal size (nodeNo*nodeNo(*epochNo))!']);
end
if size(realConn, 1) ~= size(nullConn, 1)
    error('Input args "realConn" and "nullConN" have different size on first dimension!');
end
if ~isnumeric(gammaValues) || ~isvector(gammaValues)
    error('Input arg "gammaValues" should a numeric vector!');
end
if ~isnumeric(omegaValues) || ~isvector(omegaValues)
    error('Input arg "omegaValues" should a numeric vector!');
end
if ~isnumeric(rep) || numel(rep) ~= 1 || ~ismember(rep, 10:1:10^4)
    error('Input arg "rep" should be a value in range 10:1:10^4!');
end

% user message
disp([char(10), 'Called multiCommDetectWrapper function with input args: ', ...
    char(10), 'Connectivity data set sized ', num2str(size(realConn)), ...
    char(10), 'Null connectivity data sized ', num2str(size(nullConn)), ...
    char(10), 'Range of gamma values: ', num2str(length(gammaValues)), ...
    ' values between ', num2str(min(gammaValues)), ' (min) and ', num2str(max(gammaValues)), ' (max)', ...
    char(10), 'Range of omega values: ', num2str(length(omegaValues)), ...
    ' values between ', num2str(min(omegaValues)), ' (min) and ', num2str(max(omegaValues)), ' (max)', ...
    char(10), 'No. of louvain runs for each (gamma, omega): ', num2str(rep)]);


%% Preallocate results arrays

% helpful shorthands
gNo = length(gammaValues);
oNo = length(omegaValues);
[nodeNo, ~, epochNo] = size(realConn);

% no. of communities
commNoMean = zeros(gNo, oNo);
commNoSD = zeros(gNo, oNo);
% Q
qMean = zeros(gNo, oNo);
qSD = zeros(gNo, oNo);
% partition distance (z-rand, rand, add rand)
zrandMean = zeros(gNo, oNo);
zrandSD = zeros(gNo, oNo);
% persistence
persMean = zeros(gNo, oNo);
persSD = zeros(gNo, oNo);
% consensus similarity - warning if result would be over 1GB, error if over
% 4GB
if (nodeNo*epochNo*gNo*oNo*8) >= 10^9 && (nodeNo*epochNo*gNo*oNo*8) < 4*10^9
    warning('Output consSim will be over 1GB (but below 4GB)!');
elseif (nodeNo*epochNo*gNo*oNo*8) >= 4*10^9
    error('Output consSim would be over 4GB - might be a mistake? Erroring out as a precaution!');
end
consSim = zeros(gNo, oNo, nodeNo*epochNo);

% measure elapsed time per (gamma, omega) loop
paramPairLoopTime = zeros(gNo, oNo);

% user message
disp([char(10), 'Preallocated results arrays....'])


%% Loop over parameter values

% user message
disp([char(10), 'Starting loops over gamma and omega values....'])

% clock for measuring overall elapsed time
totalClock = tic;

parfor gIdx = 1:gNo 
    g = gammaValues(gIdx);
    
    for oIdx = 1:oNo
        o = omegaValues(oIdx);
        
        % clock for measuring per-loop time
        loopClock = tic;
        
        %% Louvain!
        
        % construct multilayer connectivity matrix from realConn and
        % nullConn using getMultiLayerConnMatrix
        multiLayerConn = getMultiLayerConnMatrix(realConn, nullConn, g, o);
        
        % temp variables holding genlouvain results from all repetitions
        repS = zeros(nodeNo*epochNo, rep);
        repQ = zeros(rep, 1);
        
        % call iterated genlouvain "rep" times
        for r = 1:rep
            
            % set no specific limit; suppress output; do not force randord;
            % select random quality-increasing moves with weights given by
            % change in Q; 
            [repS(:, r), repQ(r)] = iterated_genlouvain(multiLayerConn, [], 0, [], 'moverandw');
            
        end
        
        
        %% Get stats of interest
        
        % no. of communities
        tmp = zeros(rep, 1);
        for j = 1:rep
            tmp(j) = length(unique(repS(:, j)));
        end
        commNoMean(gIdx, oIdx) = mean(tmp);
        commNoSD(gIdx, oIdx) = std(tmp);
        
        % Q
        qMean(gIdx, oIdx) = mean(repQ);
        qSD(gIdx, oIdx) = std(repQ);
        
        % consensus similarity
        [consSim(gIdx, oIdx, :), ~, pairwiseSim] = consensus_similarity(repS');
        
        % get partition distances from pairwise similarities calculated for
        % consensus similarity
        tmp = triu(pairwiseSim, 1);  % pairwiseSim is a full symmetric matrix, get the unique pairwise values = upper triangle above diagonal
        tmp = tmp(:);  % vectorize
        tmp(tmp==0)=[];  % delete the values corresponding to lower triangle and diagonal
        zrandMean(gIdx, oIdx) = mean(tmp);
        zrandSD(gIdx, oIdx) = std(tmp);

        % persistence
        persTmp = zeros(rep, 1);
        for j = 1:rep
            % reshape partition
            tmp = repS(:, j); 
            tmp = reshape(tmp, [nodeNo, epochNo]); 
            persTmp(j) = sum(sum(tmp(:,2:end) == tmp(:, 1:end-1)));
        end
        persMean(gIdx, oIdx) = mean(persTmp);
        persSD(gIdx, oIdx) = std(persTmp);
        
        % elapsed time for (gamma, omega) loop
        paramPairLoopTime(gIdx, oIdx) = toc(loopClock);
        
        
    end  % for o
    
    % user message
    disp([char(10), 'Finished all omega loops with gamma ', num2str(g)]);
    
end  % for g
       
% overall elapsed time
totalTime = toc(totalClock);

% user message
disp([char(10), 'Total time elapsed: ', num2str(totalTime)]);


%% Save and return

saveFile = 'louvainRes_alpha_stim1.mat';
save(saveFile, 'commNoMean', 'commNoSD', 'qMean', 'qSD', 'consSim', ...
    'zrandMean', 'zrandSD', 'persMean', 'persSD', 'paramPairLoopTime', ...
    'gammaValues', 'omegaValues', 'rep');


return
















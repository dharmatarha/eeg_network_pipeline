function [res, paramPairLoopTime] = multiCommDetectNullWrapper(realConn, nullConn, gammaValues, omegaValues, varargin)
%% Wrapper for repeated multilayer community detection on a randomized connectivity data set
%
% USAGE: [res, paramPairLoopTime] = multiCommDetectNullWrapper(realConn, 
%                                           nullConn, 
%                                           gammaValues, 
%                                           omegaValues,
%                                           randConn='connectionalConstrained',
%                                           rep=100, 
%                                           outputFile=[],
%                                           rawOutput= 'rawOutputNo',
%                                           method='iterated', 
%                                           randmove='moverandw', 
%                                           postprocess=[])
%
% The function explores a range of resolution parameters (spatial and
% temporal, that is, gamma and omega, respectively) for a randomized (null)
% multilayer (usually temporal) connectivity data set. The number of 
% repetitions, the genlouivan method details and output handling can all be
% controlled with input args. The function by default returns a struct with 
% a range of measures describing the partitions for each (gamma, omega) pairing.
% It  also returns the consensus similarity partitions (one partition per 
% (gamma, omega) pair).
%
% Randomization is performed inside the repetitions loop, that is, for each
% run with a give (gamma, omega) pair a newly randomized connectivity
% tensor is used.
%
% The function requires the GenLouvain toolbox, functions 
% consensus_similarity.m and zrand.m from the Network Community Toolbox, 
% our connectivity randomization functions, plus the getMultiLayerConnMatrix.m 
% and calcMultiModMatrix.m functions.
%
% PARFOR!!! Check your parallel computing settings before calling the
% function, as it runs with default parfor settings.
%
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
%               within-node inter-layer connections, that is, all 
%               inter-layer edges have the same weight. 
%
% Optional inputs:
% randConn    - Char array, one of {'connectional', 
%               'connectionalConstrained', 'nodal', 'temporal'}. Controls
%               the type of connectivity randomization we apply to obtain a
%               null model. Defaults to 'connectionalConstrained'.
% rep         - Numeric value, number of repetitions for running
%               the Louvain algorithm with given (gamma, omega) pair. 
%               Must be a positive integer in the range [10:10^4]. Defaults
%               to 100.
% outputFile  - Char array, path of a file to save the outputs into.
%               Defaults to [], that is, no output file.
% rawOutput   - Char array, one of {'rawOutputYes', 'rawOutputNo'}.
%               Controls if the exact partitions from all repetitions are
%               included in the output. There is a built-in hard upper cap 
%               of 2 GB though: if the partitions are estimated to require 
%               >2 GB memory, they are not kept internally and are not 
%               returned in the output struct. In practice, returning all 
%               partitions can only be done form small networks and few 
%               repetitions.Defaults to 'rawOutputNo'. 
% method      - Char array, one of {'iterated', 'single'}. Controls if
%               a single run of the Louvain algorithm is called or
%               the iterated version (genlouvain.m or
%               iterated_genlouvain.m from the GenLouvain toolbox). 
% randmove    - Char array, one of {'move', 'moverand', 'moverandw'}. Gets
%               passed on to genlouvain.m/iterated_genlouvain.m as randmove
%               input arg. Defaults to 'moverandw'.
% postprocess - Char array, one of {'postprocess_ordinal_multilayer',
%               'postprocess_categorical_multilayer'}. Postprocessing 
%               method for iterative_genlouvain.m, gets passed on to that 
%               function as a function handle. If left empty, there is no 
%               postprocessing. Defaults to [] (no postprocessing). Set to 
%               [] if called with method=='single'.
%
%
% Output:
% res         - Struct with the following fields:
%    commNoMean - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Mean number of communities for all (gamma, omega) pairs 
%           (averaged over repetitions).  
%    commNoSD   - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Standard deviation of number of communities for given 
%           (gamma, omega) pair. 
%    qMean      - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Mean quality for all (gamma, omega) pairs.
%    qSD      - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Standard deviation of quality for all (gamma, omega) pairs.
%    zrandMean  - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Mean partition similarity (z-Rand score) for all (gamma, omega) pairs.
%    zrandSD    - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Standard deviation of partition similarity (z-Rand score) for 
%           all (gamma, omega) pairs.
%    persMean   - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Mean persistence for all (gamma, omega) pairs. Persistence is
%           defined as the number of occasions when a node (vertex) does
%           not change its community membership from one layer to the next.
%    persSD     - Numeric matrix sized (numel(gammaValues), numel(omegaValues)). 
%           Standard deviation of persistence for all (gamma, omega) pairs.
%    consSim    - Uint16 numeric tensor sized 
%           (numel(gammaValues), numel(omegaValues), node no. X epoch no.).
%           Contains the consensus similarity partition for all (gamma,
%           omega) pairs. Warning is given if it would exceed 1 GB, not
%           calculated or included if it would exceed 2 GB.
%    part       - Uint16 numeric tensor sized
%           (numel(gammaValues), numel(omegaValues), node no. X epoch no., rep).   
%           Contains all the partitions from all repetitions for all
%           (gamma, omega) pairs. Warning is given if it would exceed 1 GB, 
%           not calculated or included if it would exceed 2 GB. In 
%           practice, it is only returned for small networks and few
%           repetitions.
% paramPairLoopTime     - Numeric matrix, sized (numel(gammaValues),
%                         numel(omegaValues)). Each value is the time the
%                         function took to calculate results for given
%                         (gamma, omega) pair.
%
% Notes:
% (1) Notice that res.consSim (and res.part if requested) are uint16, so
% can only store positive integers up to 2^16-1. This is fine for
% partitions for all of our use cases - if it does not fit your application
% surely other things would break down as well...
% (2) Fields "consSim" and "part" of the output variable might reach the
% size of 2-2 GB in certain cases. Keep that in mind. 
%

%% Input checks

% overall number of args
if ~ismembertol(nargin, 4:11)
    error(['Function multiCommDetectNullWrapper requires input args "realConN", '... 
        '"nullConn", "gammaValues" and "omegaValues" while input args "randConn", '...
        '"rep", "outputFile", rawOutput", "method", "randmove" and "postprocess" ',...
        'are optional!']);
end
% check mandatory args
if ~isnumeric(realConn) || length(size(realConn)) ~= 3 || size(realConn, 1) ~= size(realConn, 2)
    error(['Input arg "realConN" should be a 3D numeric array with the first ',...
        'two dimensions having euqal size (nodeNo*nodeNo*epochNo)!']);
end
if ~isnumeric(nullConn) || ~ismembertol(length(size(nullConn)), 2:3)  || size(nullConn, 1) ~= size(nullConn, 2)
    error(['Input arg "nullConN" should be a 2D or 3D numeric array with the first ',...
        'two dimensions having euqal size (nodeNo*nodeNo(*epochNo))!']);
end
if size(realConn, 1) ~= size(nullConn, 1)
    error('Input args "realConn" and "nullConN" have different sizes on first dimension!');
end
if length(size(nullConn)) == 3 && ~isequal(size(realConn), size(nullConn))
    error('Input arg "nulConn" is 3D but its size does not equal the size of "realConn"!');
end
if ~isnumeric(gammaValues) || ~isvector(gammaValues)
    error('Input arg "gammaValues" should be a numeric vector!');
end
if ~isnumeric(omegaValues) || ~isvector(omegaValues)
    error('Input arg "omegaValues" should be a numeric vector!');
end
% check optional args
if ~isempty(varargin)
    for v = 1:length(varargin)
        % rep
        if isnumeric(varargin{v}) && numel(varargin{v})==1 && ismembertol(varargin{v}, 10:1:10^4) && ~exist('rep', 'var')
            rep = varargin{v};
        % randConn
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'connectional', 'connectionalConstrained', 'nodal', 'temporal'}) && ~exist('randConn', 'var')
            randConn = varargin{v};            
        % rawOutput
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'rawOutputYes', 'rawOutputNo'}) && ~exist('rawOutput', 'var')
            rawOutput = varargin{v};
        % method    
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'iterated', 'single'}) && ~exist('method', 'var') 
            method = varargin{v};
        % randmove   
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'move', 'moverand', 'moverandw'}) && ~exist('randmove', 'var') 
            randmove = varargin{v};
        % postprocess    
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'postprocess_ordinal_multilayer', 'postprocess_categorical_multilayer'}) && ~exist('postprocess', 'var') 
            postprocess = varargin{v};
        % outputFile
        elseif ischar(varargin{v}) && ~exist('outputFile', 'var')
            outputFile = varargin{v};
        else
            error('At least one input arg could not be mapped to any optional input arg nicely!');
        end
    end
end
% assign default values
if ~exist('rep', 'var')
    rep = 100;
end
if ~exist('randConn', 'var')
    randConn = 'connectionalConstrained';
end
if ~exist('rawOutput', 'var')
    rawOutput = 'rawOutputNo';
end
if ~exist('method', 'var')
    method = 'iterated';
end
if ~exist('randmove', 'var')
    randmove = 'moverandw';
end
if ~exist('postprocess', 'var')
    postprocess = [];
end
if ~exist('outputFile', 'var')
    outputFile = [];    
else  % try saving a variable to the given file path
    tmp = 0; 
    try
        save(outputFile, 'tmp');
        delete(outputFile);
    catch ME
        disp([char(10), 'Failed when tested outputFile: ', outputFile,...
            ' as a valid save target!']);
        rethrow(ME);
    end
end

% further checks / transformations
% set postprocess to empty if method=='single', otherwise turn it into
% function handle
if strcmp(method, 'single')
    postprocess = [];
elseif strcmp(method, 'iterative')
    if strcmp(postprocess, 'postprocess_ordinal_multilayer')
        postprocess = @postprocess_ordinal_multilayer;
    elseif strcmp(postprocess, 'postprocess_categorical_multilayer')
        postprocess = @postprocess_categorical_multilayer;
    end
end
% turn rawOutput into boolean var
if strcmp(rawOutput, 'rawOutputYes')
    rawOutput = true;
elseif strcmp(rawOutput, 'rawOutputNo')
    rawOutput = false;
end

% user message
disp([char(10), 'Called multiCommDetectNullWrapper function with input args: ', ...
    char(10), 'Connectivity data set sized ', num2str(size(realConn)), ...
    char(10), 'Null connectivity data sized ', num2str(size(nullConn)), ...
    char(10), 'Range of gamma values: ', num2str(length(gammaValues)), ...
    ' values between ', num2str(min(gammaValues)), ' (min) and ', num2str(max(gammaValues)), ' (max)', ...
    char(10), 'Range of omega values: ', num2str(length(omegaValues)), ...
    ' values between ', num2str(min(omegaValues)), ' (min) and ', num2str(max(omegaValues)), ' (max)', ...
    char(10), 'No. of louvain runs for each (gamma, omega): ', num2str(rep),...
    char(10), 'Output file: ', outputFile,...
    char(10), 'Raw partitions in output: ', num2str(rawOutput),...
    char(10), 'Louvain method (single vs iterative): ', method,...
    char(10), 'Node selection method for Louvain algorithm: ', randmove,...
    char(10), 'Postprocessing included (valid only if iterative Louvain is requested): ', postprocess,...
    char(10)]);


%% Check the size of potentially large outputs (res.consSim & res.part)

% helpful shorthands
gNo = length(gammaValues);
oNo = length(omegaValues);
[nodeNo, ~, epochNo] = size(realConn);

% flag for calculating (and storing) consensus similarity partitions
consSimFlag = false;
% check size of consensus similarity partitions. Set flag to true if 
% size <= 2 GB, warning if >1GB. Assume uint16 (2 bytes per value)
if (nodeNo*epochNo*gNo*oNo*2) <= 2*10^9
    consSimFlag = true;
    if (nodeNo*epochNo*gNo*oNo*2) >= 10^9
        warning([char(10), 'Output res.consSim will be over 1 GB (but below 2 GB)!']);
    else
        disp([char(10), 'Output res.consSim remains under 1 GB.']);
    end
else
    warning([char(10), 'Output res.consSim would exceed 2 GB, will not be calculated and returned!']);
end

% check size of all partitions only if rawOutput was requested
% assume uint16 (2 bytes per value)
if rawOutput
    if (nodeNo*epochNo*gNo*oNo*rep*2) <= 2*10^9
        if (nodeNo*epochNo*gNo*oNo*rep*2) >= 10^9
            warning([char(10), 'Output res.part will be over 1 GB (but below 2 GB)!']);
        else
            disp([char(10), 'Output res.part remains under 1 GB.']);
        end
    else
        warning([char(10), 'Output res.part would exceed 2 GB, will not be calculated and returned!']);
        rawOutput = false;
    end
end


%% Preallocate results arrays

% no. of communities
commNoMean = zeros(gNo, oNo);
commNoSD = zeros(gNo, oNo);
% Q
qMean = zeros(gNo, oNo);
qSD = zeros(gNo, oNo);
% partition distance (z-rand)
zrandMean = zeros(gNo, oNo);
zrandSD = zeros(gNo, oNo);
% persistence
persMean = zeros(gNo, oNo);
persSD = zeros(gNo, oNo);
% consensus similarity - only if passed the size check
if consSimFlag
    consSim = zeros(gNo, oNo, nodeNo*epochNo, 'uint16');
end
% raw partition data - only if requested and passed the size check
if rawOutput
    part = zeros(gNo, oNo, nodeNo*epochNo, rep, 'uint16');
end

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
        
        % call genlouvain "rep" times
        for r = 1:rep
            
            
            % network randomization step, depending on randConn arg
            switch randConn
                case 'connectional'
                    
                case 'connectionalConstrained'
                    
                case 'nodal'
                
                case 'temporal'
                    
            end
            
            % construct multilayer connectivity matrix from realConn and
            % nullConn using getMultiLayerConnMatrix
            multiLayerConn = getMultiLayerConnMatrix(realConn, nullConn, g, o); 
            
            
            % GenLouvain specs:
            % set no specific limit; suppress output; do not force randord;
            % select quality-increasing moves according to "randmove"; set
            % no initial partition (S0); set postprocessing according to
            % "postprocess"            
            
            % select iterated or single-run version
            if strcmp(method, 'iterative')
                [repS(:, r), repQ(r)] = iterated_genlouvain(multiLayerConn, [], 0, [], randmove, [], postprocess);
            elseif strcmp(method, 'single')     
                [repS(:, r), repQ(r)] = genlouvain(multiLayerConn, [], 0, [], randmove);
            end
            
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
        
        % consensus similarity - call consensus_similarity.m with the first 
        % (consensus partition) output only if consSimFlag
        if consSimFlag
            [consSimTmp, ~, pairwiseSim] = consensus_similarity(repS');
            consSim(gIdx, oIdx, :) = uint16(consSimTmp);  % partitions are positive integers 
        else
            [~, ~, pairwiseSim] = consensus_similarity(repS');
        end
        
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
        
        % raw partition data is stored if requested and passed the memory test
        if rawOutput
            part(gIdx, oIdx, :, :) = uint16(repS);  % partitions are positive integers 
        end
        
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


%% Assign output arrays to struct fields 

% output var
res = struct;

% assign fields
res.commNoMean = commNoMean;
res.commNoSD = commNoSD;
res.qMean = qMean;
res.qSD = qSD;
res.zrandMean = zrandMean;
res.zrandSD = zrandSD;
res.persMean = persMean;
res.persSD = persSD;
% consSim and part are only assigned if relevant flags are true
if consSimFlag
    res.consSim = consSim;
end
if rawOutput
    res.part = part;
end


%% Save if file path was provided, return

if ~isempty(outputFile)
    save(outputFile, 'res', 'paramPairLoopTime', ...
       'gammaValues', 'omegaValues', 'rep', 'rawOutput',... 
       'method', 'randmove', 'postprocess');
end


return
















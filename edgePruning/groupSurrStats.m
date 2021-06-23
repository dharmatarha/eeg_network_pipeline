function [p, d] = groupSurrStats(realV, surrMu, surrSigma, varargin)
%% Statistical testing of a group of values against surrogate distributions
%
% USAGE: [p, d] = groupSurrStats(realV, 
%                               surrMu, 
%                               surrSigma, 
%                               truncateBounds=[],
%                               verbosity='verbose')
%
% The function derives the significance (p) values for the mean of a set of 
% values when tested against a set of surrogate normal distributions 
% characterized with their means (surrMu) and standard deviations 
% (surrSigma).
%
% Our frequent use case is connectivity values: a set of connectivity
% values are tested against surrogate normal distributions derived from
% phase-scrambled or otherwise randomized data.
%
% NOTE! Only the upper triangles are taken into account when supplied with
% 3D arrays containing layers of connectivity matrices.
%
% Inputs:
% realV         - Numeric array, either a vector or a 3D array. If a
%               vector, it simply contains the "real" values of some
%               measure (e.g. connectivity values). If 3D array, we expect
%               it as a multilayer connectivity tensor with dimensions 
%               [nodes X nodes X layers]. In this case, only the upper
%               triangle of each [nodes X nodes] matrix is taken into 
%               account, and values are grouped together along the 3rd
%               dimension.
% surrMu        - Numeric array, either a vector or a 3D array. Similarly 
%               to "realV", if a vector, it simply contains the mu params 
%               for the surrogate distributions of some measure (e.g. 
%               connectivity values). If 3D array, we expect it as mu 
%               params corresponding to a multilayer connectivity tensor 
%               with dimensions [nodes X nodes X layers]. In this case, 
%               only the upper triangle of each [nodes X nodes] matrix is 
%               taken into  account, and values are grouped together along
%               the 3rd dimension. Its size must match the size of "realV".
% surrSigma     - Numeric array, either a vector or a 3D array. Similarly 
%               to "surrMu", if a vector, it simply contains the sigma 
%               params for the surrogate distributions of some measure  
%               (e.g. connectivity values). If 3D array, we expect it as 
%               sigma params corresponding to a multilayer connectivity 
%               tensor with dimensions [nodes X nodes X layers]. In this case, 
%               only the upper triangle of each [nodes X nodes] matrix is 
%               taken into  account, and values are grouped together along
%               the 3rd dimension. Its size must match the size of "realV".
%
% Optional inputs:
% truncateBounds    - Numeric vector with length of 2. If supplied, the
%               normal distributions are truncated with the bounds in
%               "truncateBounds" before evaluation. The two values
%               correspond to the lower and upper bounds. Defaults to no
%               truncation (empty vector).
% verbosity         - Char array, one of {'silent', 'verbose'}. Controls if
%               the function should give user feedbacks to the command
%               line. Defaults to 'silent'.
%
% Outputs:
% p             - Numeric value or matrix, depending on the size of the
%               input arrays (vectors or 3D arrays). Contains the
%               significance values from comparing the real measurement
%               values to the surrogate normal distributions. Corrected for
%               a two-tailed test.
% d             - Numeric value or matrix, depending on the size of the
%               input arrays (vectors or 3D arrays). Contains the
%               direction of the difference between the mean real value and
%               the mean of the combined surrogate normals.
%


%% Input checks

% check no. of args
if ~ismember(nargin, 3:5)
    error('Function groupSurrStats requires input args "realV", "surrMu" and "surrSigma", while input args "truncateBounds" and "verbosity" are optional!');
end
% check mandatory args
if ~isnumeric(realV) || ~(isvector(realV) || numel(size(realV))==3)
    error('Input arg "realV" should be a numeric vector or a 3D numeric array!');
end
if ~isnumeric(surrMu) || ~(isvector(surrMu) || numel(size(surrMu))==3)
    error('Input arg "surrMu" should be a numeric vector or a 3D numeric array!');
end
if ~isnumeric(surrSigma) || ~(isvector(surrSigma) || numel(size(surrSigma))==3)
    error('Input arg "surrSigma" should be a numeric vector or a 3D numeric array!');
end
% check optional arg
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && numel(varargin{v})==2 && ~exist('truncateBounds', 'var')
            truncateBounds = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'silent', 'verbose'}) && ~exist('verbosity', 'var')
            verbosity = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to optional args "truncateBounds" and "verbosity"!');
        end
    end
end
% assign defaults
if ~exist('truncateBounds', 'var')
    truncateBounds = [];
end
if ~exist('verbosity', 'var')
    verbosity = 'verbose';
end

% sanity check, size of mandatory inputs should match
if ~isequal(size(realV), size(surrMu)) || ~isequal(size(realV), size(surrSigma))
    error('Input args "realV", "surrMu" and "surrSigma" should have the same size!');
end

% set flags
if isvector(realV)
    vectorFlag = true;
else
    vectorFlag = false;
end
if strcmp(verbosity, 'verbose')
    verbosity = true;
else
    verbosity = false;
end

% user message
if verbosity
    disp([char(10), 'Called function groupSurrStats with inputs: ',...
        char(10), 'Real (actual) values: array with dims ', num2str(size(realV)),...
        char(10), 'Surrogate mu params: array with dims ', num2str(size(surrMu)),...
        char(10), 'Surrogate  sigma params: array with dims ', num2str(size(surrSigma)),...
        char(10), 'Truncation bounds (if any): ', num2str(truncateBounds),...
        char(10), 'Verbosity: ', num2str(verbosity)]);
end


%% Case: inputs are vectors

if vectorFlag
    
    % calculate the params of the combined normal distribution
    combMu = mean(surrMu);
    combSigma = sqrt(sum((surrSigma.^2))/(numel(surrSigma)^2));
    % get significance of averaged real values
    if isempty(truncateBounds)
        p = 2 * normcdf(mean(realV), combMu, combSigma);  % two-tailed test
    else
        pd = makedist('normal', 'mu', combMu, 'sigma', combSigma);  % returns a probability distribution object
        pd = truncate(pd, truncateBounds(1), truncateBounds(2));
        p = 2 * pd.cdf(mean(realV));  % two-tailed test  
    end
    % get direction of difference
    d = (mean(realV) > combMu) * 2 - 1;  % returns 1 or -1 depending on whether the mean value is larger than the combined Mu
    
end  % if vectorFlag


%% Case: inputs are 3D arrays

if ~vectorFlag
    
    % get the number of nodes and layers in input arrays
    [nodeNo1, nodeNo2, layerNo] = size(realV);
    if nodeNo1 ~= nodeNo2
        error('The first two dimension of the input arrays should be equal!');
    end
    nodeNo = nodeNo1;

    % preallocate results vars
    p = nan(nodeNo);
    d = p;

    for node1 = 1:nodeNo 
        for node2 = 1:nodeNo   
            if node1 < node2          

                % calculate the params of the combined normal distribution
                combMu = mean(surrMu(node1, node2, :));
                combSigma = sqrt(sum((surrSigma(node1, node2, :).^2))/(layerNo^2));
                % get significance of averaged real values
                if isempty(truncateBounds)
                    p(node1, node2) = 2 * normcdf(mean(realV(node1, node2, :)), combMu, combSigma);  % two-tailed test
                else
                    pd = makedist('normal', 'mu', combMu, 'sigma', combSigma);  % returns a probability distribution object
                    pd = truncate(pd, truncateBounds(1), truncateBounds(2));
                    p(node1, node2) = 2 * pd.cdf(mean(realV(node1, node2, :)));  % two-tailed test  
                end
                % get direction of difference
                d(node1, node2) = (mean(realV(node1, node2, :)) > combMu) * 2 - 1;  % returns 1 or -1 depending on whether the mean value is larger than the combined Mu

            end  % if node1 < node2    
        end  % for node2
    end  % for node1

end  % ~if vectorFlag

return


function res = orthAmpCorr(epochData, varargin)
%% Pairwise-orthogonalized envelope correlation
%
% USAGE: res = orthAmpCorr(epochData, lpFilter=[], verbosity=1)
%
% Calculates envelope correlations on pairwise-orthogonalized ROI/channel
% data. Expects real-valued data as input in ROIs/channels X samples
% format. 
%
% As pairwise-orthogonalization is noncommutative, it is performed in both
% directions and the average of the two connectivity values is returned in
% the output arg.
%
% Mandataory input:
% epochData     - Numeric matrix where each row is a separate channel / ROI 
%               time series and columns correspond to samples. 
%
% Optional inputs:
% lpFilter      - Digital filter object as returned by e.g. designfilt
%               (part of Signal Processing Toolbox!). Used for lowpass
%               filtering the amplitude envelopes before correlations.
%               Applied via the Signal Processing Toolbox (!) version of 
%               the filter command. Defaults to [], meaning no filtering. 
% verbosity     - Verbosity. If 1, it prints to command window, 0 means
%               silence. Default is 1.
%
% Output:
% res           - Numeric matrix of connectivity values, sized
%               ROIs/channels X ROIs/channels, so that entry i,j is the
%               envelope correlation of ROIs/channels i and j. Only the
%               upper triangle above the main diagonal is populated with 
%               values, rest is NaN.
%
%
% Relevant papers:
% Coquelet et al., 2020. Comparing MEG and high-density EEG for intrinsic
% functional connectivity mapping. NeuroImage.
%


%% Input checks

if ~ismember(nargin, 1:3)
    error('Function orthAmpCorr requires input arg "epochData" while args "lpFilter" and "verbosity" are optional!');
end
if ~isnumeric(epochData) || ~ismatrix(epochData)
    error('Input arg "epochData" should be  numeric matrix!');
end
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isa(varargin{v}, 'digitalFilter') && ~exist('lpFilter', 'var')
            lpFilter = varargin{v};
        elseif ismember(varargin{v} ,[0 1]) && ~exist('verbosity', 'var')
            verbosity = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "lpFilter" or "verbosity"!');
        end
    end
end
if ~exist('lpFilter', 'var')
    lpFilter = [];
end
if ~exist('verbosity', 'var')
    verbosity = 1;
end

% check input data size
[roiNo, sampleNo] = size(epochData);
if roiNo >= sampleNo
    warning('Input data seems to have many ROIs/channels relative to the number of samples. We proceed but it is suspicious.');
end

% user message
if verbosity
    disp([char(10), 'Called orthAmpCorr on data with ', num2str(roiNo),... 
        ' ROIs/channels, each with ', num2str(sampleNo), ' samples']);
end


%% Loops through all ROI/channel pairings

epochData = epochData';

% temp var for connectivity res
res = nan(roiNo);

% vector norms (Euclidean) of all channels 
roiNorms = sqrt(sum(epochData.^2, 1));  % roiNorms is a row vector

% get envelopes for all channels
epochDataEnv = envelope(epochData);  % keep dims ROIs/channels X samples

% filter all channels if lowpass filter was provided
if ~isempty(lpFilter)
    epochDataEnv = (filter(lpFilter, epochDataEnv));
end

% loops trough ROIs
for roiIdx = 1:roiNo 
    
    % select one ROI/channel to correlate with all others, replicate it
    % roiNo times
    roiData = repmat(epochData(:, roiIdx), [1, 62]);  

    % scalar projections of data1 on all ROIs/channels
    projScalars = dot(roiData, epochData)./roiNorms;  % projScalar is 1 X roiNo
    
    % projection vectors of data1 on all ROIs/channels
    projVectors = (epochData./roiNorms).*projScalars;  
    
    % orthogonalized data1 with regards to all other channels
    orthVectors = roiData - projVectors;
    
    % get envelopes
    roiDataEnv = envelope(orthVectors);  
    
    % filter envelopes
    if ~isempty(lpFilter)
        roiDataEnv = (filter(lpFilter, roiDataEnv));
    end
    
    % get correlations
    rhos = corr(roiDataEnv, epochDataEnv);
    
    % take only the diagonal for the results matrix
    res(roiIdx, :) = diag(rhos);
    
end

% average the upper and lower triangles into the upper
% triangle
res = (triu(res, 1) + tril(res, -1)')/2;
% lower triangle is set to NaN
res(tril(true(roiNo))) = NaN;  

% user message
if verbosity
    disp(['Calculated orthAmpCorr for ', num2str(roiNo*(roiNo-1)/2), ' ROI/channel pairings']);
end


return








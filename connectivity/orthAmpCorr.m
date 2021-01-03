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

% temp var for connectivity res
res = nan(roiNo);

% loops trough ROIs
for roi1 = 1:roiNo 
    for roi2 = 1:roiNo
        if roi1 ~= roi2

            % orthogonalize dataOne with respect to dataTwo
            data1 = squeeze(epochData(roi1, :));  
            data2 = squeeze(epochData(roi2, :));
            projScalar = dot(data1, data2)/norm(data2);  % scalar projection of dataOne on dataTwo
            projVector = projScalar*(data2./norm(data2));  % projection vector of dataOne on dataTwo
            orthVector = data1 - projVector;  % orthogonalized dataOne
            % get envelope of both data 
            data1Env = envelope(orthVector);
            data2Env = envelope(data2);
            % filter envelopes
            if ~isempty(lpFilter)
                data1Env = filter(lpFilter, data1Env);
                data2Env = filter(lpFilter, data2Env);
            end
            % get amplitude correlation
            res(roi1, roi2) = corr(data1Env', data2Env');

        end  % if
    end  % for roiTwo
end  % for roiOne

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








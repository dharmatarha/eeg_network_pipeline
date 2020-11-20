function eeglabFirFilter(subFolder, subIds, firParams)
%% Function for filtering into frequency bands using EEGLAB firws filters
%
% USAGE: eeglabFirFilter(subFolder, subIds)
%
% The function defines a set of FIR filters using EEGLAB tools for
% filtering into common freuqency bands (delta, theta, alpha, beta, and 
% gamma) then applies those to subject-level EEGLAB-style data structures.
% Input args "subFolder" and "subIds" define the EEG data structs to filter
% and their path. 
%
% We expect that subject-level files are named "subjectID.mat" and contain
% a var named "EEG" which is the EEGLAB-style data struct. 
%
% FIR filters are designed on the basis of supplied params in arg
% "firParams", however, we have sensible default values for common
% frequency bands. See details below.
%
% Filtered versions of the data are saved out into "subFolder", named
% "subjectID_frequencyBand.mat". 
%
% ONLY EVER TESTED WITH KAISER WINDOWS
%
% Mandatory inputs:
% subFolder         - Char array, defines the path to the folder containing
%                   subject-level data files (*.mat files). 
% subIds            - Either a cell array with each cell containing a 
%                   subject identifier (char array) or a char array, one of
%                   {'l3_n28', 'lh_n22', 'lp1_n23', 'lp2_n25', 'lp3_n25',... 
%                   'ls_n26', 'lvid_n25'}. These char arrays refer to
%                   specific datasets with corresponding sets of subject 
%                   identifiers.
%
% Optional input:
% firParams         - Struct, with a field for each target freuqency band.
%                   E.g. firParams."alpha" contain filter params for alpha
%                   frequency band. Frequency band names must be one of
%                   {'delta', 'theta', 'alpha', 'beta', 'gamma'}.
%                   For each frequency band, the following
%                   parameters need to be defined as subfields: 
%                   (1) Cutoff frequencies as field ".cutoffs": numeric
%                   vector with two values (Hz), the first for lowpass, the
%                   second for highpass. Nan value means no filtering -
%                   e.g. [4, NaN] for delta means only lowpass filtering,
%                   no highpass.
%                   (2) Transition bandwidth as field ".transBW": numeric
%                   vector with two values (Hz), the first for lowpass, the
%                   second for highpass. 
%                   (3) Passband ripple as field ".ripple": numeric
%                   vector with two values, the first for lowpass, the
%                   second for highpass. 
%                   (4) Window type as field ".wtype": char array, one of
%                   {'kaiser', 'hamming'}. Same type is applied for low-
%                   and highpass.
%
% 
% Defaults for firParams:
%
% % DELTA
% firParams.delta.cutoffs = [4.5, NaN];
% firParams.delta.transBW = [2, NaN];
% firParams.delta.ripple = [0.002, 0.002];
% firParams.delta.wtype = 'kaiser';
% % THETA
% firParams.theta.cutoffs = [8.5, 3.5];
% firParams.theta.transBW = [2, 2];
% firParams.theta.ripple = [0.002, 0.002];
% firParams.theta.wtype = 'kaiser';
% % ALPHA
% firParams.alpha.cutoffs = [12.5, 7.5];
% firParams.alpha.transBW = [2, 2];
% firParams.alpha.ripple = [0.002, 0.002];
% firParams.alpha.wtype = 'kaiser';
% % BETA
% firParams.beta.cutoffs = [31, 11];
% firParams.beta.transBW = [3, 3];
% firParams.beta.ripple = [0.002, 0.002];
% firParams.beta.wtype = 'kaiser';
% % GAMMA
% firParams.gamma.cutoffs = [NaN, 29];
% firParams.gamma.transBW = [NaN, 6];
% firParams.gamma.ripple = [0.002, 0.002];
% firParams.gamma.wtype = 'kaiser'; 
%
%
% NOTES:
% Basic info: 
%   https://sccn.ucsd.edu/wiki/Firfilt_FAQ
%   https://github.s3.amazonaws.com/downloads/widmann/firfilt/firfilt.pdf
%


%% Input checks

% a few usefule vars
freqs = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
freqFields = {'cutoffs', 'transBW', 'ripple', 'wtype'};
definedSets = {'l3_n28', 'lh_n22', 'lp1_n23', 'lp2_n25', 'lp3_n25', 'ls_n26', 'lvid_n25', 'lv_n28'};

% check no.of args
if ~ismember(nargin, 2:3)
    error(['Function eeglabFirFilter requires input args "subFolder" and ',...
        '"subIds" while input arg "firParams" is optional!']);
end
% check mandatory args
if ~ischar(subFolder) || ~exist(subFolder, 'dir')
    error('Input arg "subFolder" is not a valid path to a folder!');
end
if ~(ischar(subIds) || iscell(subIds)) 
   error(['Input arg "subIds" must be either a cell array ',...
       'containing subject identifiers (char arrays) or a ',...
       'char array referring to a predefined dataset!']);
elseif iscell(subIds)
    if ~all(cellfun(@ischar, subIds))
        error('If input arg "subIds" is a cell array, it must contain char arrays (subject identifiers)!');
    end
elseif ischar(subIds)
    if ~ismember(subIds, definedSets)
        error(['If input arg "subIds" is a char array, it must be one of ',...
            '{"l3_n28", "lh_n22", "lp1_n23", "lp2_n25", "lp3_n25", "ls_n26", "lvid_n25", "lv_n28"}!']);
    end
end
% check optional arg
if nargin == 3
    if ~isstruct(firParams) || isempty(fieldnames(firParams)) || ~all(ismember(fieldnames(firParams), freqs))
        error(['Optional input arg "firParams" must be a struct with field names ',...
            'from {"delta", "theta", "alpha", "beta", "gamma"}!']);
    else
        fnames = fieldnames(firParams);
        for i = 1:length(fnames)
            if isempty(fieldnames(firParams.(fnames{i}))) || ~all(ismember(fieldnames(firParams.(fnames{i})), freqFields))
                error(['Input arg firParams does not have the required fields for firParams.', fnames{i}, '!']);
            end
        end
    end
else
    firParams = struct;
    for f = freqs
        firParams.(f{:}) = struct;
    end
    % DELTA
    firParams.delta.cutoffs = [4.5, NaN];
    firParams.delta.transBW = [2, NaN];
    firParams.delta.ripple = [0.002, 0.002];
    firParams.delta.wtype = 'kaiser';
    % THETA
    firParams.theta.cutoffs = [8.5, 3.5];
    firParams.theta.transBW = [2, 2];
    firParams.theta.ripple = [0.002, 0.002];
    firParams.theta.wtype = 'kaiser';
    % ALPHA
    firParams.alpha.cutoffs = [12.5, 7.5];
    firParams.alpha.transBW = [2, 2];
    firParams.alpha.ripple = [0.002, 0.002];
    firParams.alpha.wtype = 'kaiser';
    % BETA
    firParams.beta.cutoffs = [31, 11];
    firParams.beta.transBW = [3, 3];
    firParams.beta.ripple = [0.002, 0.002];
    firParams.beta.wtype = 'kaiser';
    % GAMMA
    firParams.gamma.cutoffs = [NaN, 29];
    firParams.gamma.transBW = [NaN, 6];
    firParams.gamma.ripple = [0.002, 0.002];
    firParams.gamma.wtype = 'kaiser';  
end

% standardize subFolder format with regards to last char
if subFolder(end) == '/'
    subFolder = subFolder(1:end-1);
end

% assign subIds for datasets if necessary
if ischar(subIds)
    switch subIds
        case 'l3_n28'
            subIds = {'l3_s01', 'l3_s02', 'l3_s03', 'l3_s04', 'l3_s05',...
                'l3_s06', 'l3_s07', 'l3_s08', 'l3_s09', 'l3_s10',...
                'l3_s11', 'l3_s12', 'l3_s13', 'l3_s14', 'l3_s15',...
                'l3_s16', 'l3_s17', 'l3_s18', 'l3_s19', 'l3_s20',...
                'l3_s21', 'l3_s23', 'l3_s24', 'l3_s25', 'l3_s26',...
                'l3_s27', 'l3_s28', 'l3_s29'};
        case 'lh_n22'
            subIds = {'lh_s01', 'lh_s02', 'lh_s03', 'lh_s04', 'lh_s05',...
                'lh_s06', 'lh_s07', 'lh_s08', 'lh_s09', 'lh_s10',...
                'lh_s11', 'lh_s12', 'lh_s13', 'lh_s14', 'lh_s15',...
                'lh_s16', 'lh_s17', 'lh_s18', 'lh_s19', 'lh_s20',...
                'lh_s21', 'lh_s22'};
        case 'lp1_n23'
            subIds = {'lp1_s01', 'lp1_s02', 'lp1_s03', 'lp1_s04', 'lp1_s05',...
                'lp1_s07', 'lp1_s08', 'lp1_s09',...
                'lp1_s11', 'lp1_s12', 'lp1_s13', 'lp1_s14', 'lp1_s15',...
                'lp1_s16', 'lp1_s17', 'lp1_s18', 'lp1_s19', 'lp1_s20',...
                'lp1_s21', 'lp1_s22', 'lp1_s23', 'lp1_s24', 'lp1_s25'};
        case 'lp2_n25'    
            subIds = {'lp2_s01', 'lp2_s02', 'lp2_s03', 'lp2_s04', 'lp2_s05',...
                'lp2_s06', 'lp2_s07', 'lp2_s08', 'lp2_s09', 'lp2_s10',...
                'lp2_s11', 'lp2_s12', 'lp2_s13', 'lp2_s14', 'lp2_s15',...
                'lp2_s16', 'lp2_s17', 'lp2_s18', 'lp2_s19', 'lp2_s20',...
                'lp2_s21', 'lp2_s22', 'lp2_s23', 'lp2_s24', 'lp2_s25'};    
        case 'lp3_n25'    
            subIds = {'lp3_s01', 'lp3_s02', 'lp3_s03', 'lp3_s04', 'lp3_s05',...
                'lp3_s06', 'lp3_s07', 'lp3_s08', 'lp3_s09', 'lp3_s10',...
                'lp3_s11', 'lp3_s12', 'lp3_s13', 'lp3_s14', 'lp3_s15',...
                'lp3_s16', 'lp3_s17', 'lp3_s18', 'lp3_s19', 'lp3_s20',...
                'lp3_s21', 'lp3_s22', 'lp3_s23', 'lp3_s24', 'lp3_s25'};   
        case 'ls_n26'
            subIds = {'ls_s01', 'ls_s02', 'ls_s03', 'ls_s04', 'ls_s05',...
                'ls_s06', 'ls_s07', 'ls_s08', 'ls_s09', 'ls_s10',...
                'ls_s11', 'ls_s12', 'ls_s13', 'ls_s14', 'ls_s15',...
                'ls_s16', 'ls_s17', 'ls_s18', 'ls_s19', 'ls_s20',...
                'ls_s21', 'ls_s22', 'ls_s23', 'ls_s24', 'ls_s25',...
                'ls_s26'};    
        case 'lvid_n25'
            subIds = {'lvid_s01', 'lvid_s02', 'lvid_s03', 'vid_s04', 'lvid_s05',...
                'lvid_s06', 'lvid_s07', 'lvid_s08', 'lvid_s09', 'lvid_s10',...
                'lvid_s11', 'lvid_s12', 'lvid_s13', 'lvid_s14', 'lvid_s15',...
                'lvid_s16', 'lvid_s17', 'lvid_s18', 'lvid_s19', 'lvid_s20',...
                'lvid_s21', 'lvid_s22', 'lvid_s23', 'lvid_s24', 'lvid_s25'};    
        case 'lv_n28'
            subIds = {'lv_s01', 'lv_s02', 'lv_s03', 'lv_s04', 'lv_s05',...
                'lv_s06', 'lv_s07', 'lv_s08', 'lv_s09', 'lv_s10',...
                'lv_s11', 'lv_s12', 'lv_s13', 'lv_s14', 'lv_s15',...
                'lv_s16', 'lv_s17', 'lv_s18', 'lv_s19', 'lv_s20',...
                'lv_s21', 'lv_s22', 'lv_s23', 'lv_s24', 'lv_s25',...
                'lv_s26', 'lv_s27', 'lv_s28'};     
    end  % switch subIds
end  % if ischar(subIds)

% user feedback
disp([newline, 'Called eeglabFirFilter with input args: ',...
    newline, 'Subject data folder: ', subFolder,...
    newline, 'Subject identifiers: ']);
disp(subIds);
disp('FIR filter parameters: ');
freqNames = fieldnames(firParams);
for i = 1:length(freqNames)
    disp(firParams.(freqNames{i}));
end


%% Get sampling rate from first file

firstFile = [subFolder, '/', subIds{1}, '.mat'];
tmp = load(firstFile);
Fs = tmp.EEG.srate;

% user feedback
disp([newline, 'Sampling rate of EEG data is ', num2str(Fs), ' in the first data file at ', firstFile, '.',...
    newline, 'We assume that all EEG data has the same sampling rate!']);


%% EEGLAB init

% user feedback
disp([newline, 'EEGLAB init,,,']);

% try to init eeglab functions
try 
    eeglab nogui;
catch ME
    warning('Could not init / load EEGLAB functions. Calling "eeglab nogui" resulted in the following error:');
    rethrow(ME);
end


%% Define filters

% get cell array for frequency bands of interest
freqNames = fieldnames(firParams);

% user feedback
disp([newline, 'Preparing filters for the following frequency bands:']);
disp(freqNames);

% loop over frequency bands
for i = 1:length(freqNames)
    currentFreq = freqNames{i};
    
    % get normalized cutoffs
    firParams.(currentFreq).cutoffsNormed = firParams.(currentFreq).cutoffs./(Fs/2);
    
    %% lowpass
    
    % if lowpass filtering is requested
    if ~isnan(firParams.(currentFreq).cutoffs(1))
        
        % define filter params based on existing fields
        if strcmp(firParams.(currentFreq).wtype, 'kaiser')
            firParams.(currentFreq).lpBeta = kaiserbeta(firParams.(currentFreq).ripple(1));  % get beta for kaiser window
            firParams.(currentFreq).lpOrder = firwsord(firParams.(currentFreq).wtype, Fs, firParams.(currentFreq).transBW(1), firParams.(currentFreq).ripple(1));  % estimate filter order
            firParams.(currentFreq).lpWindow = windows(firParams.(currentFreq).wtype, firParams.(currentFreq).lpOrder+1, firParams.(currentFreq).lpBeta);  % get window
        elseif strcmp(firParams.(currentFreq).wtype, 'hamming')
            [firParams.(currentFreq).lpOrder, firParams.(currentFreq).ripple(1)]= firwsord(firParams.(currentFreq).wtype, Fs, firParams.(currentFreq).transBW(1));  % estimate filter order
            firParams.(currentFreq).lpWindow = windows(firParams.(currentFreq).wtype, firParams.(currentFreq).lpOrder+1);  % get window
        end
        firParams.(currentFreq).lpCoeffs = firws(firParams.(currentFreq).lpOrder, firParams.(currentFreq).cutoffsNormed(1), firParams.(currentFreq).lpWindow);  % lowpass filter coeffs
    end

    
    %% highpass

    % if highpass filtering is requested
    if ~isnan(firParams.(currentFreq).cutoffs(2))     
        
        % define filter params based on existing fields
        if strcmp(firParams.(currentFreq).wtype, 'kaiser')
            firParams.(currentFreq).hpBeta = kaiserbeta(firParams.(currentFreq).ripple(2));  % get beta for kaiser window
            firParams.(currentFreq).hpOrder = firwsord(firParams.(currentFreq).wtype, Fs, firParams.(currentFreq).transBW(2), firParams.(currentFreq).ripple(2));  % estimate filter order
            firParams.(currentFreq).hpWindow = windows(firParams.(currentFreq).wtype, firParams.(currentFreq).hpOrder+1, firParams.(currentFreq).hpBeta);  % get window
        elseif strcmp(firParams.(currentFreq).wtype, 'hamming')
            [firParams.(currentFreq).hpOrder, firParams.(currentFreq).ripple(2)]= firwsord(firParams.(currentFreq).wtype, Fs, firParams.(currentFreq).transBW(2));  % estimate filter order
            firParams.(currentFreq).hpWindow = windows(firParams.(currentFreq).wtype, firParams.(currentFreq).hpOrder+1);  % get window
        end
        firParams.(currentFreq).hpCoeffs = firws(firParams.(currentFreq).hpOrder, firParams.(currentFreq).cutoffsNormed(2), 'high', firParams.(currentFreq).hpWindow);  % highpass filter coeffs
    end  
    
    
    %% Convolve filters if both are requested
    
    if all(~isnan(firParams.(currentFreq).cutoffs))
        if isequal(firParams.(currentFreq).lpOrder, firParams.(currentFreq).hpOrder)
            firParams.(currentFreq).lphpCoeffs = conv(firParams.(currentFreq).lpCoeffs, firParams.(currentFreq).hpCoeffs, 'same');
        else
            error(['It is not fully straightforward how to convolve the low- and highpass filters for band ',... 
                currentFreq, ' due to different orders. ',...
                'Choose the same transition bandwidth values and ripple values for both filters']);
        end
    end
    
    % user feedback
    disp([newline, 'Filter params for ', currentFreq, ':']);
    disp(firParams.(currentFreq));
    
    
end  % for i = 1:length(freqNames)
    
% user feedback
disp([newline, 'Prepared all filters']);


%% Apply filters to all data

% check for / create subfolders for frequency bands
% loop over frequency bands
for i = 1:length(freqNames)
    currentFreq = freqNames{i};
    % check for existing folder for current freq
    if ~exist([subFolder, '/', currentFreq], 'dir')
        mkdir([subFolder, '/', currentFreq]);
    end    
end

% loop over data files
for s = 1:length(subIds)
    
    subjectClock = tic;
    
    % user feedback
    disp([newline, 'Filtering data from subject ', subIds{s}, '...']);
    
    % load subject-level data
    subFile = [subFolder, '/', subIds{s}, '.mat'];
    data = load(subFile);
    subEEG = data.EEG;

    % sanity checks
    if ~isequal(subEEG.srate, Fs) || isempty(subEEG.data)
        error(['Bad sampling rate or missing data at ', subFile, '!']);
    end

    % treat a problem with earlier processing of epoch-data - missing
    % nbchan field
    if ~isfield(subEEG, 'nbchan')
        subEEG.nbchan = size(subEEG.data, 1);
    end
    
    % user feedback
    disp(['Loaded data for subject ', subIds{s}, ', data size is:']);
    disp(size(subEEG.data));
    
    % loop over frequency bands
    for i = 1:length(freqNames)
        currentFreq = freqNames{i};
        
        % user feedback
        disp(['Filtering for ', currentFreq, '...']);
        
        % copy subject data for filtering
        tmpData = subEEG;
        
        % save file for given freq
        freqSaveP = [subFolder, '/', currentFreq, '/', subIds{s}, '_', currentFreq, '.mat'];
        
        % convolved filter if both low- and highpass were requested
        if all(~isnan(firParams.(currentFreq).cutoffs))
            tmpData = firfilt(tmpData, firParams.(currentFreq).lphpCoeffs);
            
        % if only lowpass filter was requested
        elseif ~isnan(firParams.(currentFreq).cutoffs(1))
            tmpData = firfilt(tmpData, firParams.(currentFreq).lpCoeffs);
        
        % if only highpass filter was requested    
        elseif ~isnan(firParams.(currentFreq).cutoffs(2))
            tmpData = firfilt(tmpData, firParams.(currentFreq).hpCoeffs);
        end
        
        % save out filtered data
        tmpData.history = ['Filtered to ', currentFreq, '. See filter details in var "firParams".'];  % add note about filtering
        EEG = tmpData;
        save(freqSaveP, 'EEG', 'firParams');
        
    end  % for i

    % user feedback - report elapsed time
    disp(['Filtering subject''s data took ',... 
        num2str(round(toc(subjectClock), 2)), ' secs']);
    
end  % for s


return















function [f1, f2] = filterResCmp(data, filteredData, fs)
%% Helper function to compare power spectra of raw vs filtered data
%
% USAGE: [f1, f2] = filterResCmp(data, filteredData, fs)
%
% The function plots the power spectrum densities of both data samples.
%
% Inputs:
% data          - Numeric vector, real.
% filteredData  - Numeric vector, real. Same length as "data".
% fs            - Numeric value, sampling rate in Hz.
%
% Outputs:
% f1            - Figure handle for power spectrum density of "data"
% f2            - Figure handle for power spectrum density of "filteredData"


%% Input checks

if nargin ~= 3 
    error('Function filterResCmp requires input args "data" and "filteredData"!');
end
if ~isnumeric(data) || ~isvector(data) || ~isreal(data)
    error('Input arg "data" should be a numeric real vector!');
end
if ~isnumeric(filteredData) || ~isvector(filteredData) || ~isreal(filteredData)
    error('Input arg "data" should be a numeric real vector!');
end
if ~isequal(length(data), length(filteredData))
    error('Input args "data" and "filteredData" should be vectors with equal length!');
end
if ~isnumeric(fs) || numel(fs)~=1
    error('Input arg "fs" should be a numeric value (sampling rate in Hz)!');
end


%% Get power spectrum densities

% length in samples
L = length(data);

% get fft
fft1 = fft(data);
fft2 = fft(filteredData);
% get only half of it 
fft1_part = fft1(1:floor(L/2)+1);
fft2_part = fft2(1:floor(L/2)+1);
% get amplitude values 
fft1_amp = abs(fft1_part);
fft2_amp = abs(fft2_part);
% get power
fft1_pow = fft1_amp.^2;
fft2_pow = fft2_amp.^2;
% multiply by two for all components excluding the DC (zero) and Nquist frequency 
% in order to correct for taking only half of the fft output 
fft1_pow(2:end-1) = fft1_pow(2:end-1)*2;
fft2_pow(2:end-1) = fft2_pow(2:end-1)*2;
% normalize with time for density
psd1 = fft1_pow/(fs*L);
psd2 = fft2_pow/(fs*L);

% frequencies for the fft
freqs = 0:fs/L:fs/2;


%% Plot densities

% data
f1 = figure;
plot(freqs, 10*log10(psd1));
title('Periodogram for "data"');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
% filtered data
f2 = figure;
plot(freqs, 10*log10(psd2));
title('Periodogram for "filteredData"');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


return


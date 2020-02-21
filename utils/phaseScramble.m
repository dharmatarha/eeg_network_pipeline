function yPhaseRand = phaseScrambleMatrix(y, imagThresh)
%% Phase-scrambling function
% 
% USAGE: yPhaseRand = phaseScramble(y, imagThresh=0.01)
%
% Scrambles (randomizes) the phases of the input matrix. It works either
% with a column vector input (single variable / time series) or a matrix
% where each column is treated as a variable / time series.
%
% Logic:
% Performs fft, swaps the phases of fourier components to uniform random
% values then calculates inverse fft.
%
% Intended for real numbers, as we only deal with half of the fft spectrum.
%
% Input: 
% y             - Numerical matrix (real), input data, each column is a
%           separately phase scrambled time series
% imagThresh    - Threshold for detecting errors. Due to numerical
%           imprecision, the pipeline fft -> random phase injection -> ifft
%           often yields complex results for real input, where the
%           imaginary part is really small. Small however is a relative
%           term - this argument controls the magnitude of imaginary part 
%           at which the function returns with an error. Defaults to 0.001.
%
% Output:
% yPhaseRand    - numerical matrix (real), data with scrambled phases
%
% Notes:
% (1) An easy test of the function is given below. The phase scrambled data
% should have the same magnitude (amplitude) of DFFT components as the
% original (a negligible difference exists due to numerical inaccuracies):
%   data = randn(1000, 10);
%   dataRand = phaseScrambleMatrix(data);
%   % difference of fft amplitudes
%   fftAmpDiff = abs(fft(data))-abs(fft(dataRand));
%   % print the maximum differences for each column
%   disp(['Maximum difference between corresponding FFT amplitude components, ',...
%   'for each variable:']);
%   disp(max(fftAmpDiff));
%


%% Input checks

if nargin == 1
    imagThresh = 0.001;
end
if ~ismatrix(y) || ~isreal(y)
    error('Function phaseScramble expects a numerical matrix of reals as input!');
end

% columns are variables - give a warning if there seems to be more
% variables than time points
if size(y, 2) > size(y, 1)
    warning(['There seems to be more time series than data points in ',...
        'each series - are you sure columns are separate time series?'])
end


%% Phase-scrambling

% time series length
L = size(y, 1);
% no. of time series
nTS = size(y, 2);

% get fft
yfft = fft(y);
% get only half of it (exclude also zero component)
yfft_part = yfft(2:floor(L/2)+1, :);
% get amplitude values
yfft_amp = abs(yfft_part);

% random phases sampled randomly from uniform distr. between -pi and +pi
randomphases = (rand([floor(L/2), nTS])*2*pi)-pi;

% get cartesian coordinates with new phases and old amplitude values
[cosCoord, sinCoord] = pol2cart(randomphases, yfft_amp);

% turn coordinates into fft components of phase scrambled data
phaseRandFFT_part = complex(cosCoord, sinCoord);
if mod(L, 2) == 0 % concatenate differently depending on even/odd data
    phaseRandFFT = [yfft(1, :); phaseRandFFT_part(1:end-1, :); yfft(L/2+1, :); flipud(conj(phaseRandFFT_part(1:end-1, :)))];
else
    phaseRandFFT = [yfft(1, :); phaseRandFFT_part; flipud(conj(phaseRandFFT_part))];
end

% inverse FFT
yPhaseRand = ifft(phaseRandFFT);


%% Correct for numerical imprecision

% Numerical methods introduce a minimal imprecision, and thus do rarely 
% yield real data after the conversion. The remaining imaginary part 
% should be really small though and we can just get rid of it. But first 
% let's test the magnitude of the remaining imaginary part

if any(any(abs(imag(yPhaseRand))>imagThresh))
    error('Stg went wrong at phase randomization, the resulting data is complex!');
else
    yPhaseRand = real(yPhaseRand);
end


return
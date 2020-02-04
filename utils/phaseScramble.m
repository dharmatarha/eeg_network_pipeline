function yPhaseRand = phaseScramble(y, imagThresh)
%% Phase-scrambling function
% 
% USAGE: yPhaseRand = phaseScramble(y, imagThresh=0.01)
%
% Scrambles (randomizes) the phases of the input vector. 
% Performs fft, swaps the phases of fourier components to random values,
% then calculated inverse fft.
%
% Intended for real numbers, as we only deal with half of the fft spectrum.
%
% Input: 
% y             - Numerical vector (real), input data
% imagThresh    - Threshold for detecting errors. Due to numerical
%           imprecision, the pipeline fft -> random phase injection -> ifft
%           often yields complex results for real input, where the
%           imaginary part is really small. Small however is a relative
%           term - this argument controls the magnitude of imaginary part 
%           at which the function returns with an error. Defaults to 0.01.
%
% Output:
% yPhaseRand    - numerical vector (real), data with scrambled phases
%

%% Input checks

if nargin == 1
    imagThresh = 0.01;
end
if ~isvector(y) || ~isreal(y)
    error('Function phaseScramble expects a numberical vector of reals as input!');
end

% need column vector as input
if ~iscolumn(y)
    y = y';
end


%% Phase-scrambling

% data length
L = length(y);

% get fft
yfft = fft(y);
% get only half of it (exclude also zero component)
yfft_part = yfft(2:floor(L/2)+1);
% get amplitude values
yfft_amp = abs(yfft_part);

% random phases
randomphases = (rand([floor(L/2), 1])*2*pi)-pi;

% get cartesian coordinates with new phases and old amplitude values
[cosCoord, sinCoord] = pol2cart(randomphases, yfft_amp);

% turn coordinates into fft components of phase scrambled data
phaseRandFFT_part = complex(cosCoord, sinCoord);
if mod(L, 2) == 0 % concatenate differently depending on even/odd data
    phaseRandFFT = [yfft(1); phaseRandFFT_part; flipud(conj(phaseRandFFT_part(1:end-1)))];
else
    phaseRandFFT = [yfft(1); phaseRandFFT_part; flipud(conj(phaseRandFFT_part(1:end)))];
end

% inverse FFT
yPhaseRand = ifft(phaseRandFFT);


%% Correct for numerical imprecision

% Numerical methods introduce a minimal imprecision, and thus do rarely 
% yield real data after the conversion. The remaining imaginary part 
% should be really small though and we can just get rid of it. But first 
% let's test the magnitude of the remaining imaginary part

if any(abs(imag(yPhaseRand))>imagThresh)
    error('Stg went wrong at phase randomization, the resulting data is complex!');
else
    yPhaseRand = real(yPhaseRand);
end


return
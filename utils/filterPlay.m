% Play around with different FIR filters


%% Generate data - known frequency sinusoid + white noise

% Sampling rate
Fs = 1000;     
dt = 1/Fs;   
% Length of sample
StopTime = 5;
% Time vector
t = (0:dt:StopTime-dt)';
% Target frequency
Fc = 60;   
% Data
x = sin(2*pi*Fc*t);
% Add noise
xn = x + randn([length(t), 1]);


%% Get phase of original signal

analyticSignal = hilbert(x);
xphase = angle(analyticSignal);


%% Check spectrum of signal and noisy signal

% length of time series data
L = length(x);
% get fft
xfft = fft(x);
% get only half of it (exclude also zero component)
xfft_part = xfft(2:floor(L/2)+1, :);
% get amplitude values 
xfft_amp = abs(xfft_part);
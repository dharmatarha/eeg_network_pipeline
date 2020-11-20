% Play around with different FIR filters


%% Generate data - known frequency sinusoid + white noise

% sampling rate
Fs = 1000;     
dt = 1/Fs;   
% length of sample
stopTime = 5;  % in secs
L = Fs*stopTime;  % in samples
% time vector
t = (0:dt:stopTime-dt)';
% signal frequency
signalF = 60;   
% generate sample
x = sin(2*pi*signalF*t);
% only noise
n = randn([length(t), 1]);
% add noise to signal
xn = x + n;


%% Get phase of original signal

analyticSignal = hilbert(x);
xphase = angle(analyticSignal);


%% Check spectrum of signal, noise and noisy signal

% get fft
xfft = fft(x);
% get only half of it 
xfft_part = xfft(1:floor(L/2)+1);
% get amplitude values 
xfft_amp = abs(xfft_part);
% get power
xfft_pow = xfft_amp.^2;
% multiply by two for all components excluding the DC (zero) and Nquist frequency 
% in order to correct for taking only half of the fft output 
xfft_pow(2:end-1) = xfft_pow(2:end-1)*2;
% normalize with time for density
xpsd = xfft_pow/(Fs*L);

% frequencies for the fft
freq = 0:Fs/L:Fs/2;

% plot
figure;
plot(freq, 10*log10(xpsd));
title('Periodogram for sinusiod sample');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% same for noisy data
xnfft = fft(xn);
xnfft_pow = abs(xnfft(1:floor(L/2+1))).^2;
xnfft_pow(2:end-1) = xnfft_pow(2:end-1)*2;
xnpsd = xnfft_pow/(Fs*L);
figure;
plot(freq, 10*log10(xnpsd));
title('Periodogram for sinusiod sample with white noise');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% same for noise
nfft = fft(n);
nfft_pow = abs(nfft(1:floor(L/2+1))).^2;
nfft_pow(2:end-1) = nfft_pow(2:end-1)*2;
npsd = nfft_pow/(Fs*L);
figure;
plot(freq, 10*log10(npsd));
title('Periodogram for white noise');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


%% Get FIR filters using EEGLAB functions

% init eeglab / put eeglab functions into path
eeglab nogui;

% lowpass

% We query optimal Kaiser window parameters for given passpand ripple and
% transition bandwidth (Hz)
lpRipple = 0.004;
lpBW = 3;
lpWtype = 'kaiser';
lpCutoff = signalF + 2;
lpCutoffNorm = lpCutoff/(Fs/2);  % normalized freq to rad/sample
% get Kaiser beta
lpBeta = kaiserbeta(lpRipple);
% estimate order
lpOrder = firwsord(lpWtype, Fs, lpBW, lpRipple);
% get window
lpWindow = windows(lpWtype, lpOrder+1, lpBeta);
% get filter
lpCoeffs = firws(lpOrder, lpCutoffNorm, lpWindow);


% highpass

hpRipple = 0.004;
hpBW = 3;
hpWtype = 'kaiser';
hpCutoff = signalF - 2;
hpCutoffNorm = hpCutoff/(Fs/2);  % normalized freq to rad/sample
% get Kaiser beta
hpBeta = kaiserbeta(hpRipple);
% estimate order
hpOrder = firwsord(hpWtype, Fs, hpBW, hpRipple);
% get window
hpWindow = windows(hpWtype, hpOrder+1, hpBeta);
% get filter
hpCoeffs = firws(hpOrder, hpCutoffNorm, 'high', hpWindow);


%% Filter signal and its noisy version

% lowpass

% sinusoid signal
xeeg = struct;
xeeg.srate = Fs;
xeeg.data = x';
xeeg.trials = 1;
xeeg.pnts = L;
xeeg.event = [];

% signal + noise
xneeg = struct;
xneeg.srate = Fs;
xneeg.data = xn';
xneeg.trials = 1;
xneeg.pnts = L;
xneeg.event = [];

% noise
% signal + noise
neeg = struct;
neeg.srate = Fs;
neeg.data = n';
neeg.trials = 1;
neeg.pnts = L;
neeg.event = [];


% lowpass
% x_lp = firfilt(xeeg, lpCoeffs, length(xeeg.data));
% xn_lp = firfilt(xneeg, lpCoeffs, length(xneeg.data));
% n_lp = firfilt(neeg, lpCoeffs, length(neeg.data));
x_lp = firfilt(xeeg, lpCoeffs);
xn_lp = firfilt(xneeg, lpCoeffs);
n_lp = firfilt(neeg, lpCoeffs);

% highpass

% x_lphp = firfilt(x_lp, hpCoeffs, length(x_lp.data));
% xn_lphp = firfilt(xn_lp, hpCoeffs, length(xn_lp.data));
% n_lphp = firfilt(n_lp, hpCoeffs, length(n_lp.data));
x_lphp = firfilt(x_lp, hpCoeffs);
xn_lphp = firfilt(xn_lp, hpCoeffs);
n_lphp = firfilt(n_lp, hpCoeffs);


% filtering in one step (convolution)
lphpCoeffs = conv(lpCoeffs, hpCoeffs, 'same');
x_conv = firfilt(xeeg, lphpCoeffs);
xn_conv = firfilt(xneeg, lphpCoeffs);
n_conv = firfilt(neeg, lphpCoeffs);


%% Check filtering results

% get fft
% xfft = fft(x_lphp.data);
xfft = fft(x_conv.data);
% get only half of it 
xfft_part = xfft(1:floor(L/2)+1);
% get amplitude values 
xfft_amp = abs(xfft_part);
% get power
xfft_pow = xfft_amp.^2;
% multiply by two for all components excluding the DC (zero) and Nquist frequency 
% in order to correct for taking only half of the fft output 
xfft_pow(2:end-1) = xfft_pow(2:end-1)*2;
% normalize with time for density
xpsd = xfft_pow/(Fs*L);

% frequencies for the fft
freq = 0:Fs/L:Fs/2;

% plot
figure;
plot(freq, 10*log10(xpsd));
title('Periodogram for sinusiod sample after filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% same for noisy data
xnfft = fft(xn_lphp.data);
% xnfft = fft(xn_conv.data);
xnfft_pow = abs(xnfft(1:floor(L/2+1))).^2;
xnfft_pow(2:end-1) = xnfft_pow(2:end-1)*2;
xnpsd = xnfft_pow/(Fs*L);
figure;
plot(freq, 10*log10(xnpsd));
title('Periodogram for sinusiod sample with white noise after filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% same for noise
% nfft = fft(n_lphp.data);
nfft = fft(n_conv.data);
nfft_pow = abs(nfft(1:floor(L/2+1))).^2;
nfft_pow(2:end-1) = nfft_pow(2:end-1)*2;
npsd = nfft_pow/(Fs*L);
figure;
plot(freq, 10*log10(npsd));
title('Periodogram for white noise after filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


%% Get phase of filtered signal and noisy signal

analyticSignal = hilbert(x_lphp.data');
filt_xphase = angle(analyticSignal);

analyticSignal = hilbert(xn_lphp.data');
filt_xnphase = angle(analyticSignal);

analyticSignal = hilbert(xn_conv.data');
filt_xnconvphase = angle(analyticSignal);

analyticSignal = hilbert(n_lphp.data');
filt_nphase = angle(analyticSignal);



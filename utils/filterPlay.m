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
% add noise
xn = x + randn([length(t), 1]);


%% Get phase of original signal

analyticSignal = hilbert(x);
xphase = angle(analyticSignal);


%% Check spectrum of signal and noisy signal

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


%% Get FIR filters using EEGLAB functions

% init eeglab / put eeglab functions into path
eeglab nogui;

% lowpass

% We query optimal Kaiser window parameters for given passpand ripple and
% transition bandwidth (Hz)
lpRipple = 0.004;
lpBW = 3;
lpWtype = 'kaiser';
lpCutoff = 62;
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
hpCutoff = 58;
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
xeeg = firfilt(xeeg, lpCoeffs, length(xeeg.data));

% signal + noise
xneeg = struct;
xneeg.srate = Fs;
xneeg.data=xn';
xneeg.trials = 1;
xneeg.pnts = L;
xneeg.event = [];
xneeg = firfilt(xneeg, lpCoeffs, length(xneeg.data));


% highpass

xeeg = firfilt(xeeg, hpCoeffs, length(xeeg.data));
xneeg = firfilt(xneeg, hpCoeffs, length(xneeg.data));


%% Check filtering results

% get fft
xfft = fft(xeeg.data);
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
xnfft = fft(xneeg.data);
xnfft_pow = abs(xnfft(1:floor(L/2+1))).^2;
xnfft_pow(2:end-1) = xnfft_pow(2:end-1)*2;
xnpsd = xnfft_pow/(Fs*L);
figure;
plot(freq, 10*log10(xnpsd));
title('Periodogram for sinusiod sample with white noise after filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


%% Get phase of filtered signal and noisy signal

analyticSignal = hilbert(xeeg.data');
filt_xphase = angle(analyticSignal);

analyticSignal = hilbert(xneeg.data');
filt_xnphase = angle(analyticSignal);


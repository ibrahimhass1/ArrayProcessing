% Sample signal
Fs = 16000; % Sampling frequency
t = 0:1/Fs:1; % Time vector

% Parameters for spectrogram
windowSize = 320; % Window size for STFT
overlap = 160; % Overlap length between windows
nfft = 1024; % FFT length

% Compute and display the spectrogram
spectrogram(s_t./sum(s_t.^2), windowSize, overlap, nfft, Fs, 'yaxis');
colorbar; % Add colorbar for intensity scale

%Specify the color range
cmin = -160; % Minimum color value in dB
cmax = -80; % Maximum color value in dB
caxis([cmin cmax]);

title('Spectrogram of noisy signal (normalized)')
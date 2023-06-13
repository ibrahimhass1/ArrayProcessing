%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', ...
    1, 'speech_shaped', 0);

%% Obtained Short-Time Fourier transform of time domains signals
s = STFT(s_t, 16000, 20, 50, 'hamming');
x = STFT(x_t, 16000, 20, 50, 'hamming');


%% Apply processing
%
%
%



%% Reconstruct signal, by taking the inverse Short-Time Fourier transform
s_hat = iSTFT(s(1, :, :), 16000, 20, 50);





%% Play audio
Fs = 16000;
soundsc(s_hat(Fs:5*Fs), Fs)



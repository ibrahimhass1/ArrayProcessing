%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', ...
    0, 'white', 0);

%% Obtained Short-Time Fourier transform of time domains signals
s = STFT(s_t, 16000, 20, 50, 'hamming');
x = STFT(x_t, 16000, 20, 50, 'hamming');


%% Apply processing

% Set smoothing parameter
alpha = .9;

% Store most recent estimate of autocorrelation matrix per bin 
Rx = zeros([size(x, 2) size(x, 1) size(x, 1)]);

% Iterate over all frames
for l = 1:size(x, 3)
    l
    % Update estimate for each frequency bin
    for k = 1:size(x, 2)
        Rx(k, :, :) = alpha.*squeeze(Rx(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
        % Obtain acoestic transfer function from Rx 
        [U, ~] = eig(squeeze(Rx(k, :, :)));
        a = U(:, 1);
        w = a./(a'*a);
        s(1, k, l) = w'*x(:, k, l);
    end
end




%% Reconstruct signal, by taking the inverse Short-Time Fourier transform
s_hat = iSTFT(s(1, :, :), 16000, 20, 50);





%% Play audio
Fs = 16000;
soundsc(s_hat(Fs:5*Fs), Fs)



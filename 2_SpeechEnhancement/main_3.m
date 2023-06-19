%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', 0);

% s_t = s_t(1:100000);
% x_t = x_t(1:100000);
%% Obtained Short-Time Fourier transform of time domains signals
s = STFT(s_t, 16000, 20, 50, 'hamming');
x = STFT(x_t, 16000, 20, 50, 'hamming');

% Store clean s as reference for SNR calculations
s_clean = s;

% Add spatial white noise for each frequency bin
sigma_n = .0025;
for k=1:size(x, 2)
    x(:, k, :) = squeeze(x(:, k, :)) + sigma_n.*(randn([size(x, 1) size(x, 3)]) + 1i* randn([size(x, 1) size(x, 3)]))./sqrt(2);
end

%% Apply processing

% Set smoothing parameter
alpha = .1;

% Store most recent estimate of autocorrelation matrix per bin 
Rx = zeros([size(x, 2) size(x, 1) size(x, 1)]);
Rn = zeros([size(x, 2) size(x, 1) size(x, 1)]);


Rn_root = zeros([size(x, 2) size(x, 1) size(x, 1)]);
R_tilde = zeros([size(x, 2) size(x, 1) size(x, 1)]);
Rs = zeros([size(x, 2) size(x, 1) size(x, 1)]);

fft_coef_sum = zeros([size(x, 3) 1]);
SNR = zeros([1 size(x, 3)]);

% Iterate over all frames
for l = 1:size(x, 3)
    % Voice Activity Detection, apply thresholding on 
    fft_coef_sum(l) = sum(sum(abs(x(:, :, l)), 2).^2);
end
threshold = .3*mean(fft_coef_sum);

% Iterate over all frames
for l = 1:size(x, 3)
% for l=1:10
    l
    % Voice Activity Detection, apply thresholding the check for each frame
    % if it is either speech or noise
    fft_coef_sum(l) = sum(sum(abs(x(:, :, l)), 2).^2);

    for k = 1:size(x, 2)
        % Speech active, update Rx, keep Rn

        if fft_coef_sum(l) > threshold
            % Update Rx
            Rx(k, :, :) = alpha.*squeeze(Rx(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
            % Keep Rn
            Rx(k, :, :) = Rx(k, :, :);

            Rn_root(k, :, :) = inv(sqrtm(squeeze(Rn(k, :, :))));
            X_tilde = squeeze(Rn_root(k,:,:))*x(:, k, l);

            R_tilde(k,:,:) = X_tilde*X_tilde';

            [U, V] = eig(squeeze(R_tilde(k, :, :)));

            sigma_s = (V(1,1)-V(2,2));
            a = U(:, 1);

            Rs_tilde = a*sigma_s*a';
            Rs(k,:,:) = sqrtm(squeeze(Rn(k, :, :)))*Rs_tilde*sqrtm(squeeze(Rn(k, :, :)));
            

        % Speech inactive, update Rn, keep Rx
        else
            % Update Rn
            Rn(k, :, :) = alpha.*squeeze(Rn(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
            % Keep Rx
            Rx(k, :, :) = Rx(k, :, :);

            Rn_root(k, :, :) = inv(sqrtm(squeeze(Rn(k, :, :))));
            X_tilde = squeeze(Rn_root(k,:,:))*x(:, k, l);

            R_tilde(k,:,:) = X_tilde*X_tilde';

            [U, V] = eig(squeeze(R_tilde(k, :, :)));

            sigma_s = (V(1,1)-V(2,2));
            a = U(:, 1);

            Rs_tilde = a*sigma_s*a';
            Rs(k,:,:) = sqrtm(squeeze(Rn(k, :, :)))*Rs_tilde*sqrtm(squeeze(Rn(k, :, :)));

        [U, V] = eig(squeeze(Rs(k, :, :)));
        a = U(:, 1);

       % w = a./(a'*a); 

        %% MVDR beamformer
        Rx_inv = inv(squeeze(Rx(k, :, :)));
        w_mvdr = (Rx_inv*a)/(a'*Rx_inv*a);  

        %% Multi-channel wiener 
        %w = (Rx_inv*a)/(a'*Rx_inv*a);  
        w_wf = sigma_s /(sigma_s+inv(a'*Rx_inv*a)) * w_mvdr; 

        %% Apply beamformer
        s(1, k, l) = w_wf'*x(:, k, l);

        end

    end
    


    % Calculate clean signal power
    clean_signal = sum(abs(s_clean(1, :, l)).^2);

    % Calculate beamformer output power
    output_signal = sum(abs(s(1, :, l)).^2);

    % Calculate SNR
    SNR(l) = clean_signal/(abs(clean_signal-output_signal));
end

s_hat = iSTFT(s(1, :, :), 16000, 20, 50);

nanIndices = isnan(s_hat);
s_hat_2 = s_hat(~nanIndices);
MSE = mean(abs(s_t(1:length(s_hat_2))' - s_hat_2).^2);


figure
plot(s_hat)
title('Beamformer output')
xlabel("sample index")

figure
plot(s_t)
title("Clean target signal")
xlabel("sample index")





% 
% 
















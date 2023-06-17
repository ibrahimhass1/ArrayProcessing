%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', 0);

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


fft_coef_sum = zeros([size(x, 3) 1]);
SNR = zeros([1 size(x, 3)]);

% Iterate over all frames
for l = 1:size(x, 3)
    % Voice Activity Detection, apply thresholding on 
    fft_coef_sum(l) = sum(sum(abs(x(:, :, l)), 2).^2);
end
threshold = .1*mean(fft_coef_sum);

% Iterate over all frames
for l = 1:size(x, 3)
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
            Rn(k, :, :) = Rn(k, :, :);   
         
        % Speech inactive, update Rn, keep Rx
        else
            % Update Rn
            Rn(k, :, :) = alpha.*squeeze(Rn(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
            % Keep Rx
            Rx(k, :, :) = Rx(k, :, :);
        end
       
        % Obtain acoestic transfer function from estiamted Rx 
        [U, V] = eig(squeeze(Rx(k, :, :)));
        sigma_s = (V(1,1)-V(2,2));
        a = U(:, 1);
        Rx_inv = inv(squeeze(Rx(k, :, :)));

        
        %% Uncomment desired beamformer
        
        %% Delay-and-sum 
        w = a./(a'*a); 
        
        %% MVDR beamformer
        %w = (Rx_inv*a)/(a'*Rx_inv*a);  

        %% Multi-channel wiener 
        %w = (Rx_inv*a)/(a'*Rx_inv*a);  
        %w = sigma_s /(sigma_s+inv(a'*Rx_inv*a)) * w; 
        
        %% Apply beamformer
        s(1, k, l) = w'*x(:, k, l);
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






















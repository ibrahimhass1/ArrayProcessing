%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', ...
    0, 'white', .0002);



%% Obtained Short-Time Fourier transform of time domains signals
s = STFT(s_t, 16000, 20, 50, 'hamming');
x = STFT(x_t, 16000, 20, 50, 'hamming');
size(x)


% Store clean s as reference for SNR calculations
s_clean = s;




sigma_n = .0025;
% Add spatial white noise for each frequency bin
for k=1:size(x, 2)
    x(:, k, :) = squeeze(x(:, k, :)) + sigma_n.*(randn([size(x, 1) size(x, 3)]) + 1i* randn([size(x, 1) size(x, 3)]))./sqrt(2);
  % x(:, k, :) = squeeze(x(:, k, :)) + wgn(4, 624, 1, 'linear', 'complex');

end
% 
% s = squeeze(s);
% for k=1:size(x, 2)
%    s(k, :);
%    random = randn([1 624]);
%    s(k, :) = s(k, :) + sigma_n.*(random)./sqrt(2);
% end

%% Apply processing

% Set smoothing parameter
alpha = .5;

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
    %(l/size(x, 3))*100
    % Voice Activity Detection, apply thresholding on 
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
       
        % Determine ATF
        %Obtain acoestic transfer function from Rx 
        [U, V] = eig(squeeze(Rx(k, :, :)));
        sigma_s = (V(1,1)-V(2,2));
        a = U(:, 1);
        Rx_inv = inv(squeeze(Rx(k, :, :)));

        % Delay-and-sum beamformer
        w = a./(a'*a);  % 0.0029 MSE
        
        % MVDR beamformer
      %  w = (Rx_inv*a)/(a'*Rx_inv*a);  % 0.0030 MSE

        % 
        %w = sigma_s /(sigma_s+inv(a'*Rx_inv*a)) * w; % 0.0030 MSE

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
title('S hat')

figure
plot(x_t(1, :))
title("X_t")



% %% Reconstruct signal, by taking the inverse Short-Time Fourier transform
% s_hat = iSTFT(x(1, :, :), 16000, 20, 50);
% 
% 
% 
% 
% 
% %% Play audio
%Fs = 16000;
%soundsc(s_hat(Fs:5*Fs), Fs)
% 
% 
% 
% 
% 
% 
% 





















% as = zeros([size(x, 3) size(x, 2) 1])
% % Iterate over all frames
% for l = 1:size(x, 3)
%     l
%     % Update estimate for each frequency bin
%     for k = 1:size(x, 2)
%         Rx(k, :, :) = alpha.*squeeze(Rx(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
%         % Obtain acoestic transfer function from Rx 
%         [U, ~] = eig(squeeze(Rx(k, :, :)));
%         a = U(:, 1);
%         w = a./(a'*a);
%         s(1, k, l) = w'*x(:, k, l);
%     end
% end
% 
% 
% 
% 
% %% Reconstruct signal, by taking the inverse Short-Time Fourier transform
% s_hat = iSTFT(x(1, :, :), 16000, 20, 50);
% 
% 
% 
% 
% 
% %% Play audio
% Fs = 16000;
% soundsc(s_hat(Fs:5*Fs), Fs)



%% Obtain clean speech signal, and generate noise microphone signals
[s_t, x_t] = genData('data\clean_speech.wav', 'data\babble_noise.wav', ...
    1, 'white', .0002);



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
alpha = .1;

% Store most recent estimate of autocorrelation matrix per bin 
Rx =.01* ones([size(x, 2) size(x, 1) size(x, 1)]);
Rn = ones([size(x, 2) size(x, 1) size(x, 1)]);


fft_coef_sum = zeros([size(x, 3) 1]);

SNR = zeros([1 size(x, 3)]);

% Iterate over all frames
for l = 1:size(x, 3)
    % Voice Activity Detection, apply thresholding on 
    fft_coef_sum(l) = sum(sum(abs(x(:, :, l)), 2).^2);
end
threshold = .3*mean(fft_coef_sum);
kl = 0
pp=0
% Iterate over all frames
for l = 1:size(x, 3)
    %(l/size(x, 3))*100
    % Voice Activity Detection, apply thresholding on 
    fft_coef_sum(l) = sum(sum(abs(x(:, :, l)), 2).^2);

    l
    
    for k = 1:size(x, 2)
        % Pre-whitening
        
        x_wiggle = inv(Rn_sqrt) * x(:, k, l);



        % Speech active, update Rx, keep Rn
        if fft_coef_sum(l) > threshold
            pp = pp+1;
            % Update Rx
            Rx(k, :, :) = alpha.*squeeze(Rx(k, :, :)) + (1-alpha).* x_wiggle* x_wiggle';
            % Keep Rn
            Rn(k, :, :) = Rn(k, :, :);  
         
        % Speech inactive, update Rn, keep Rx
        else
            kl = kl+1;
            % Update Rn
            Rn(k, :, :) = alpha.*squeeze(Rn(k, :, :)) + (1-alpha).* x_wiggle* x_wiggle';
            Rn_sqrt = squeeze(Rn(k, :, :))^(1/2);
            % Keep Rx
            Rx(k, :, :) = Rx(k, :, :);
        end


        % Determine ATF


        %Obtain acoestic transfer function from Rx 
        [U, V] = eig(squeeze(Rx(k, :, :)));

        V1 = zeros(4);
        V1(1,1) = V(1,1) -1;
      %  sigma_s = (V(1,1)-V(2,2));
        U1 = zeros(4);
        U1(:,1) = U(:, 1); 
        
        [U, V] = eig(Rn_sqrt*U1*V1*U1'*Rn_sqrt);

        a = U(:, 1);


       % a = Rn_sqrt * U(:, 1);
        Rx_inv = inv(squeeze(Rx(k, :, :)));

        % Delay-and-sum beamformer
        w = a./(a'*a);  % 0.0029 MSE
        
        % MVDR beamformer
      %  w = (Rx_inv*a)/(a'*Rx_inv*a);  % 0.0030 MSE
%
        % 
    %    w = sigma_s /(sigma_s+inv(a'*Rx_inv*a)) * w; % 0.0030 MSE

        s(1, k, l) = w'*x_wiggle;
    end

    % Calculate clean signal power
%     clean_signal = sum(abs(s_clean(1, :, l)).^2);
% 
%     % Calculate beamformer output power
%     output_signal = sum(abs(s(1, :, l)).^2);
%     
%     % Calculate SNR
%     SNR(l) = clean_signal/(abs(clean_signal-output_signal));
end

s_hat = iSTFT(s(1, :, :), 16000, 20, 50);

% nanIndices = isnan(s_hat);
% s_hat_2 = s_hat(~nanIndices);

%MSE = mean(abs(s_t(1:length(s_hat_2))' - s_hat_2).^2);


figure
plot(s_hat)
title('Received signal at one of the microphones')
xlabel("Time(s)")

% figure
% plot((1:length(s_t))./Fs, s_t)
% title("Clean target signal")
% xlabel("Time(s)")

%plot(10*log10(SNR))

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

soundsc(s_hat(Fs:5*Fs), Fs)



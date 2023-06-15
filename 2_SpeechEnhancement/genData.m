% Returns target speech signal, and received microphone signals in time
% domain. 
function [s, x] = genData(speech_file, interference_file, ...
    number_of_interferers, noise_type, snr)

    %% Load data set
    % Retreive data from audio file
    [s_t, Fs] = audioread(speech_file);
    [s_i, Fs] = audioread(interference_file);

    % Make sure impulse responses are in workspace
    load('data\impulse_responses.mat')

    %% Construct received signal in time domain
    % Construct the noiseless microphone signals by convolving the speech
    % signal with the room impulse response
    for i=1:size(h_target, 1)
        x_t(i, :) = conv(s_t, h_target(i, :));
    end
    

    % Detect parts with active speech
     

    signal_mean = mean(mean(abs(x_t), 2));

    signal_power = mean(mean(x_t(abs(x_t) > signal_mean).^2, 2));

    noise_power = signal_power/(10^(snr/10));



    %% Add interference in time domain
    % Add interference by convolving the interference speech signals with
    % the impulse responses, max 4 interferers
    for i=1:size(h_target, 1)
        if number_of_interferers > 0
            x_i = conv(s_i, h_inter1(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 1
            x_i = conv(s_i(100:end), h_inter2(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 2
            x_i = conv(s_i(200:end), h_inter3(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 3
            x_i = conv(s_i(300:end), h_inter4(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
    end

    %% Add noise in time domain 
    switch lower(noise_type)
        case 'white'
            [n_t] = noise_power.*randn(4, length(x_t));
        case 'non-stationary'
            [n_t, ~] = audioread('data\aritificial_nonstat_noise.wav');
        case 'speech_shaped'
            [n_t, ~] = audioread('data\Speech_shaped_noise.wav');
        otherwise
            error('PARAMETER ERROR: noise type not reconized')
    end

%     for i=1:size(h_target, 1)
%         x_t(i, :) = x_t(i, :); % + n_t(1: length(x_t))';
%     end


    x_t = x_t + n_t;



    %% Take STFT of clean and received signal
    
    s = s_t';
    x = x_t(:, 1:length(s_t));





   

end


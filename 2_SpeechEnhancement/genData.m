function [s, x] = genData(speech_file, number_of_interferers, ...
    noise_type, snr)

    %% Load data set
    % Retreive data from audio file
    [s_t, Fs] = audioread(speech_file);
    [s_i, Fs] = audioread('data\babble_noise.wav');

    % Make sure impulse responses are in workspace
    load('data\impulse_responses.mat')

    %% Construct received signal in time domain
    % Construct the noiseless microphone signals by convolving the speech
    % signal with the room impulse response
    for i=1:size(h_target, 1)
        x_t(i, :) = conv(s_t, h_target(i, :));
    end

    %% Add interference in time domain 
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
            [n_t] = randn(length(x_t), 1);
        case 'non-stationary'
            [n_t, ~] = audioread('data\aritificial_nonstat_noise.wav');
        case 'speech_shaped'
            [n_t, ~] = audioread('data\Speech_shaped_noise.wav');
        otherwise
            error('PARAMETER ERROR: noise type not reconized')
    end

    for i=1:size(h_target, 1)
        x_t(i, :) = x_t(i, :) + n_t(1: length(x_t))';
    end


    %% Take STFT of clean and received signal



    x = x_t;
    s=s_t;





    
    %% Set processing parameters
    % frame duration in ms
    frame_duration = 20;                 
    frame_size = floor(frame_duration*Fs/1000);
    
    % overlap-and-add frame percentage shift
    overlap_percentage = 50;
    overlap = frame_size*overlap_percentage/100;
    
    % define window type
    window = hamming(frame_size);

    % Determine unscaled noise from noise type argument
%     switch lower(noise_type)
%         case 'white'
%             [n_t, ~] = randn(length(s_t));
%         case 'non-stationary'
%             [n_t, ~] = audioread('data\aritificial_nonstat_noise.wav');
%         case 'babble'
%             [n_t, ~] = audioread('data\babble_noise.wav');
%         case 'speech_shaped'
%             [n_t, ~] = audioread('data\Speech_shaped_noise.wav');
%         otherwise
%             error('PARAMETER ERROR: noise type not reconized')
%     end


     
%     for i=1:size(h_target, 1)
%         x_t(i, :) = conv(s_t, h_target(i, :));
%     end
    
    i = 1;
    
    
%     
%     % Overlap-and-add 
%     for n=1:overlap:length(s_t)-frame_size
%     s_windowed = s_t(n:n+frame_size-1).*window;
%     s(:, i) = fft(s_windowed);
%     i = i+1;
%     end


end


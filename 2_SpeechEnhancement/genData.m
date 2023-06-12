function [s, x] = genData(speech_file, number_of_interferers, ...
    noise_type, snr)


    %% Load data set
    % Retreive data from audio file
    [s_t, Fs] = audioread(speech_file);
    % Make sure impulse responses are in workspace
    load('data\impulse_responses.mat')
    
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
    switch lower(noise_type)
        case 'white'
            [n_t, ~] = randn(length(s_t));
        case 'non-stationary'
            [n_t, ~] = audioread('data\aritificial_nonstat_noise.wav');
        case 'babble'
            [n_t, ~] = audioread('data\babble_noise.wav');
        case 'speech_shaped'
            [n_t, ~] = audioread('data\Speech_shaped_noise.wav');
        otherwise
            error('PARAMETER ERROR: noise type not reconized')
    end


     
    for i=1:size(h_target, 1)
        x_t(i, :) = conv(s_t, h_target(i, :));
    end
    
    i = 1;
    
    
    
    % Overlap-and-add 
    for n=1:overlap:length(s_t)-frame_size
    s_windowed = s_t(n:n+frame_size-1).*window;
    s(:, i) = fft(s_windowed);
    i = i+1;
    end


end


% Returns target speech signal, and received microphone signals in time
% domain. 
function [s, x] = genData(speech_file, interference_file, number_of_interferers)

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

    %% Add interference in time domain
    % Add interference by convolving the interference speech signals with
    % the impulse responses, max 4 interferers
    for i=1:size(h_target, 1)
        if number_of_interferers > 0
            x_i = conv(.1*s_i, h_inter1(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 1
            x_i = conv(.1*s_i(100:end), h_inter2(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 2
            x_i = conv(.1*s_i(200:end), h_inter3(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
        
        if number_of_interferers > 3
            x_i = conv(.1*s_i(300:end), h_inter4(i, :));
            x_t(i, :) = x_t(i, :) + x_i(1:length(x_t))';
        end
    end
    % Reformat signals
    s = s_t';
    x = x_t(:, 1:length(s_t));
end


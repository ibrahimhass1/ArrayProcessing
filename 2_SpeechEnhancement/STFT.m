% frame duration in ms   
function [x] = STFT(x_t, Fs, frame_duration, ...
    overlap_percentage, window_type)

    %% Set processing parameters            
    frame_size = floor(frame_duration*Fs/1000);
    overlap = floor(frame_size*overlap_percentage/100);

    % Calculate window function
    switch lower(window_type)
        case 'bartlett'
            window = bartlett(frame_size);
        case 'gausswin'
            window = gausswin(frame_size);
        case 'hamming'
            window = hamming(frame_size);
        case 'hann'
            window = hann(frame_size);
        otherwise
            error('PARAMETER ERROR: window type not reconized')
    end

    %% Apply windowing, take FFT return as 2D array

    % zero padding??
    
    x = zeros([size(x_t, 1) frame_size floor((size(x_t, 2)-overlap)/overlap)]);

    if size(x_t, 1) > 1
        for j=1:size(x_t, 1)
            i = 1;
            % Overlap-and-add 
            for n=1:overlap:length(x_t)-frame_size
                x_windowed = x_t(j, n:n+frame_size-1).*window';
                x(j, :, i) = fft(x_windowed);
                i = i+1;
            end
        end
    else
        i = 1;
        % Overlap-and-add 
        for n=1:overlap:length(x_t)-frame_size
            x_windowed = x_t(n:n+frame_size-1).*window';
            x(1, :, i) = fft(x_windowed);
            i = i+1;
        end
    end
end


% 
%          for j=1:size(x_t, 1)
%         i = 1;
%         % Overlap-and-add 
%         for n=1:overlap:length(x_t)-frame_size
%             x_windowed = x_t(n:n+frame_size-1).*window;
%             x(:, i) = fft(x_windowed);
%             i = i+1;
%         end

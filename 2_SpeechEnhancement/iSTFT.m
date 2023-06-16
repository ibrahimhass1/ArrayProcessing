% frame duration in ms   
function [x_t] = iSTFT(x, Fs, frame_duration, overlap_percentage)

    %% Set processing parameters            
    window_size = floor(frame_duration*Fs/1000);
    overlap = floor(window_size*overlap_percentage/100);

   

    %% Apply windowing, take FFT return as 2D array
    
    x_t = zeros([overlap*(size(x, 3)+1), 1]);
    x = squeeze(x);
    % Overlap-and-add 
    for l = 1:size(x, 2)
        x_t((l-1)*overlap+1: (l-1)*overlap+window_size) = real(x_t((l-1)*overlap+1: (l-1)*overlap+window_size) + ifft(x(:,l)));
    end
end

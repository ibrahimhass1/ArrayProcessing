function [MVDR_estimates] = MVDR_estimates(x,a)

    % Set smoothing parameter
    alpha = .9;
    
    % Store most recent estimate of autocorrelation matrix per bin 
    Rx = zeros([size(x, 2) size(x, 1) size(x, 1)]);
    
    % Iterate over all frames
    for l = 1:size(x, 3)
        % Update estimate for each frequency bin
        for k = 1:size(x, 2)
            Rx(k, :, :) = alpha.*squeeze(Rx(k, :, :)) + (1-alpha).*x(:, k, l)*x(:, k, l)';
            % Obtain acoestic transfer function from Rx 
        end
    end
    
    
    
    MVDR_estimates = inv(Rx).*a/(a'.*Rx^(-1).*a);

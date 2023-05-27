function Xsm = temporal_smoothing(X,m)
    
    % Check if smoothing factor is valid, that is a positive non-zero
    % integer
    if or((m <= 0), not(mod(m, 1) == 0))
      error('Invalid smoothing factor');
    end

    % Smooth X by a factor m, that is calculate moving average over m samples
    Xsm = zeros(size(X));
  
    
    M = size(X, 1);
    N = size(X, 2);
    % Zero padd
    X = [zeros([size(X, 1) m-1]), X];
    
    % Apply smooting, store result in Xsm
    for j=1:M
        for i=1:N
            Xsm(j, i) = mean(X(j, i:i+m-1));
        end
    end
end



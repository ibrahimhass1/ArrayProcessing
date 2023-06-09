function s_hat = zeroForcingReceiver(h, X, N)

% Determine oversampling factor
P = length(h)

% Construct H
nul = zeros([P 1]);
H = [h, nul; nul, h];

% Calculate estimate for S
S_hat = pinv(H)*X;

% Determine symbol estimates, first and last we only have
% one estimate available
s_hat = zeros([N 1]);
s_hat(1) = S_hat(1,1);
% As for the rest average over the 2 estimates
s_hat(N) = S_hat(2, end);
s_hat(2:N-1) = (S_hat(1, 2:end) + S_hat(2, 1:end-1))./2;

% Translate to nearest qpsk symbol
for k = 1:N
    if real(s_hat(k)) > 0 && imag(s_hat(k)) > 0
        s_hat(k) = (1 + 1i) / sqrt(2);
    elseif real(s_hat(k)) <= 0 && imag(s_hat(k)) > 0
        s_hat(k) = (-1 + 1i) / sqrt(2);
    elseif real(s_hat(k)) > 0 && imag(s_hat(k)) <= 0
        s_hat(k) = (1 - 1i) / sqrt(2);
    else
        s_hat(k) = (-1 - 1i) / sqrt(2);
    end
end







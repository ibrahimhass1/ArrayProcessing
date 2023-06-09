function [S_hat_wiener, s_hat, s_raw] = WienerReceiver(X,h,sigma, N)

nul = zeros([length(h),1]);
H = [h,nul;nul,h];

[m, n] = size(H*H');
fprintf('The matrix has %d rows and %d columns.\n', m, n);
% +  wgn(2*length(h),2*length(h),0.5^2)
Wf = pinv(H*H'+ sigma^2*eye(2*length(h)))*H;

S_hat_wiener = Wf'*X;

S_hat = S_hat_wiener;


% Determine symbol estimates, first and last we only have
% one estimate available
s_hat = zeros([N 1]);
s_hat(1) = S_hat(1,1);
% As for the rest average over the 2 estimates
s_hat(N) = S_hat(2, end);
s_hat(2:N-1) = (S_hat(1, 2:end) + S_hat(2, 1:end-1))./2;

% Also return raw estimates
s_raw = s_hat;

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
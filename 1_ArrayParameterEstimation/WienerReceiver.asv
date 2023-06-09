function S_hat_wiener = WienerReceiver(X,h)

nul = zeros([length(h),1]);
H = [h,nul;nul,h];

[m, n] = size(H*H');
fprintf('The matrix has %d rows and %d columns.\n', m, n);
% +  wgn(2*length(h),2*length(h),0.5^2)
Wf = pinv(H*H'+  wgn(2*length(h),2*length(h),0.5^2) )*H;

S_hat_wiener = Wf'*X;
    
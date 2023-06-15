%% Q3, implement ESPRIT algorithm
% Generate data using the genData function
[M, N, Delta, theta, f, SNR] = deal(3, 20, 0.5, [-20; 30], [0.1; 0.12], 60);
[X, A, S] = genData(M, N, Delta, theta, f, SNR);
ESPRIT_freq(X, 2)
%% Zero-forcing beamformer with A is known

% % Obtain estimates of source directions, sort such that signal estimate
% % will have the same order
% DoA = sort(esprit(X,2));
% % Construct zero-forcing beamformer
% phi = 2*pi*Delta*sin(deg2rad(DoA));
% A_hat = exp(1i*(0:M-1)'.*phi');
% W = pinv(A_hat);
% % Apply beamformer to obtain an estimate of S
% S_hat = W*X;
% 
% %% Zero-forcing beamformer with S is known
% 
% % Obtain estimates of source directions, sort such that signal estimate
% % will have the same order
% freq = sort(espritfreq(X,2));
% % Construct zero-forcing beamformer
% S_hat =  exp(1i*2*pi*freq*(0:N-1));
% % Apply beamformer to obtain an estimate of S
% A_hat = X*pinv(S_hat);
% 
% i=1
% for theta=-90:90
%     exponent = 2*pi*Delta*sin(deg2rad(theta));
%     a = exp(1i*(0:M-1)'.*exponent);
%     y(i, :) = abs(W*a);
%     i = i+1;
% end
% 
% figure
% hold on
% plot(-90:90, y(:, 1))
% plot(-90:90, y(:, 2))
% xlabel("Theta(degree)")
% ylabel("|y(theta)|")
% hold off
% 
% N = 500;
% s = ones([N 1]);
% %h = [1; 1; -1; -1; 1; 1; -1; -1];
% h = [1; -1; 1; -1];
% qpsk_symbols = [-1-1i, 1-1i, -1+1i, 1+1i]./sqrt(2);
% 
% for i=1:N
%     s(i) = qpsk_symbols(randi([1 4]));
% end
% 
% P = 4;
% sigma = 4.5;
% X = gendata_conv(s,P,N,sigma);
% %x = gendata_conv_2(s,P,N,sigma);
% rank(X)
% 
% 
% [S, s_hat, s_raw] = zeroForcingReceiver(h, X, N);
% [S_wf, s_hat_wf, s_raw_wf] = WienerReceiver(X, h,sigma, N);
% figure;
% hold on;
% plot(s_raw_wf, '.', 'color', 'blue', 'MarkerSize', 5);
% plot(qpsk_symbols, 'x', 'color', 'red', 'MarkerSize', 15, 'LineWidth', 3);
% title("Raw estimates of Wiener receiver (N=500, P=8)")
% xlabel("Re(s)")
% ylabel("Im(s)")
% xlim([-1.5,  1.5]);
% ylim([-1.5,  1.5]);
% hold off;
% 
% mse_wiener = mean(abs(s-s_raw_wf).^2)
% 
% figure;
% hold on;
% plot(s_raw, '.', 'color', 'blue', 'MarkerSize', 5);
% plot(qpsk_symbols, 'x', 'color', 'red', 'MarkerSize', 15, 'LineWidth', 3);
% title("Raw estimates of Zero-Forcing receiver (N=500, P=8)")
% xlabel("Re(s)")
% ylabel("Im(s)")
% xlim([-1.5,  1.5]);
% ylim([-1.5,  1.5]);
% hold off;
% 
% mse_zero = mean(abs(s-s_raw).^2)
% 


% 
% 
% %% Wiener receiver
% 
% R_x = X*X';
% r = X*x;
% 
% w = inv(R_x);




% [X, A, S] = genDataKR(M, N, Delta, theta, f, SNR);
% size(X)
% % Call the esprit function
% AoA = esprit(X, 2)
% freq = espritfreq(X, 2)

% Apply temporal smoothing, with a factor m
% X = eye(4)
% m = 1;
% Xs = temporal_smoothing(X,m)


% % Q2, effect of number of samples on the singular values
% sv_N20 = 0;
% sv_N40 = 0;
% sv_N60 = 0;
% sv_N100 = 0;
% 
% for i=1:100
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_N20=sv_N20+svd(X);
%     [X,A,S] = genData(5,40,.5,[-20; 30],[.1;.3],20);
%     sv_N40=sv_N40+svd(X);
%     [X,A,S] = genData(5,80,.5,[-20; 30],[.1;.3],20);
%     sv_N60=sv_N60+svd(X);
%     [X,A,S] = genData(5,160,.5,[-20; 30],[.1;.3],20);
%     sv_N100=sv_N100+svd(X);
% end
% sv_N20 = sv_N20/100
% sv_N40 = sv_N40/100
% sv_N60 = sv_N60/100
% sv_N100 = sv_N100/100
% 
% figure
% hold on
% plot((sv_N100), '.', 'MarkerSize', 20)
% plot((sv_N60), '.', 'MarkerSize', 20)
% plot((sv_N40), '.', 'MarkerSize', 20)
% 
% plot((sv_N20), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Singular value")
% legend('N=160', 'N=80', 'N=40', 'N=20')
% title("Singular values of X for different number of samples")
% hold off
% 
% 10^mean(log10(sv_N100)-log10(sv_N20))
% 10^mean(log10(sv_N60)-log10(sv_N20))
% 10^mean(log10(sv_N40)-log10(sv_N20))
% 
% % Q2, effect of number of antennas on the singular values
% sv_A15 = 0;
% sv_A10 = 0;
% sv_A5 = 0;
% 
% for i=1:100
%     [X,A,S] = genData(15,20,.5,[-20; 30],[.1;.3],20);
%     sv_A15 =sv_A15 +svd(X);
%     [X,A,S] = genData(10,20,.5,[-20; 30],[.1;.3],20);
%     sv_A10 =sv_A10 +svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_A5=sv_A5+svd(X);
% end
% sv_A15 = sv_A15/100
% sv_A10 = sv_A10/100
% sv_A5 = sv_A5/100
% % sv_A1 = sv_A1/100
% 
% 
% figure
% hold on
% plot((sv_A15), '.', 'MarkerSize', 20)
% plot((sv_A10), '.', 'MarkerSize', 20)
% plot((sv_A5), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:15) 
% ylabel("Singular value")
% legend('M=15', 'M=10', 'M=5')
% title("Singular values of X for different number of antennas")
% hold off
% 
% 10^mean(log10(sv_A15(1:5))-log10(sv_A10(1:5)))
% 10^mean(log10(sv_A10(1:5))-log10(sv_A5(1:5)))
% 10^mean(log10(sv_A15(1:5))-log10(sv_A5(1:5)))

%% Q2, effect of the angle between the sources on the singular values
% sv_A50 = 0;
% sv_A10 = 0;
% sv_A1 = 0;
% 
% for i=1:100
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_A50 =sv_A50 +svd(X);
%     [X,A,S] = genData(5,20,.5,[20; 30],[.1;.3],20);
%     sv_A10=sv_A10+svd(X);
%     [X,A,S] = genData(5,20,.5,[29; 30],[.1;.3],20);
%     sv_A1=sv_A1+svd(X);
% end
% sv_A50 = sv_A50/100
% sv_A10 = sv_A10/100
% sv_A1 = sv_A1/100
% 
% 
% figure
% hold on
% plot((sv_A50), '.', 'MarkerSize', 20)
% plot((sv_A10), '.', 'MarkerSize', 20)
% plot((sv_A1), '.', 'MarkerSize', 20)
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Singular value")
% legend('50 deg', '10 deg', '1 deg')
% title("Singular values of X for different angles between sources")
% hold off


%% Q2, effect of the source frequency difference on the singular values;
% sv_f1 = 0;
% sv_f05 = 0;
% sv_f025 = 0;
% sv_f01 = 0;
% sv_f001 = 0;
% 
% for i=1:100
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.2;.3],20);
%     sv_f1=sv_f1+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.25;.3],20);
%     sv_f05=sv_f05+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_f025 =sv_f025 +svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.29;.3],20);
%     sv_f01=sv_f01+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.299;.3],20);
%     sv_f001=sv_f001+svd(X);
% end
% sv_f1 = sv_f1/100
% sv_f05 = sv_f05/100
% sv_f025 = sv_f025/100
% sv_f01 = sv_f01/100
% sv_f001 = sv_f001/100
% 
% 
% figure
% hold on
% plot((sv_f1), '.', 'MarkerSize', 20)
% plot((sv_f05), '.', 'MarkerSize', 20)
% plot((sv_f025), '.', 'MarkerSize', 20)
% plot((sv_f01), '.', 'MarkerSize', 20)
% plot((sv_f001), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Singular value")
% legend('.1', '.05', '0.25', '.01', '.001')
% title("Singular values of X for different frequency differences")
% hold off





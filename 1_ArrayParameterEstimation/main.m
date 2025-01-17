%% Signal model, Q2 - Effect number of samples on singular values
% sv_N20 = 0;
% sv_N40 = 0;
% % Average over 100 realizations
% for i=1:100
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_N20=sv_N20+svd(X);
%     [X,A,S] = genData(5,40,.5,[-20; 30],[.1;.3],20);
%     sv_N40=sv_N40+svd(X);
% end
% sv_N20 = sv_N20/100;
% sv_N40 = sv_N40/100;
% 
% figure
% hold on
% plot((sv_N40), '.', 'MarkerSize', 20)
% plot((sv_N20), '.', 'MarkerSize', 20)
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Singular value")
% legend('N=40', 'N=20')
% title("Singular values of X for different number of samples")
% hold off

%% Signal model, Q2 - Effect number of antennas on singular values
% sv_A10 = 0;
% sv_A5 = 0;
% % Average over 100 realizations
% for i=1:100
%     [X,A,S] = genData(10,20,.5,[-20; 30],[.1;.3],20);
%     sv_A10 =sv_A10 +svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_A5=sv_A5+svd(X);
% end
% sv_A10 = sv_A10/100;
% sv_A5 = sv_A5/100;
% figure
% hold on
% plot((sv_A10), '.', 'MarkerSize', 20)
% plot((sv_A5), '.', 'MarkerSize', 20)
% xlabel("Singular value index")
% xticks(1:15) 
% ylabel("Singular value")
% legend('M=10', 'M=5')
% title("Singular values of X for different number of antennas")
% hold off

%% Signal model, Q2 - Effect AoA difference on singular values
% sv_A50 = 0;
% sv_A10 = 0;
% sv_A1 = 0;
% % Average over 100 realizations
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

%% Signal model, Q2 - Effect source frequency difference on singular values
% sv_f1 = 0;
% sv_f05 = 0;
% sv_f025 = 0;
% sv_f01 = 0;
% sv_f001 = 0;
% % Average over 100 realizations
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
% sv_f1 = sv_f1/100;
% sv_f05 = sv_f05/100;
% sv_f025 = sv_f025/100;
% sv_f01 = sv_f01/100;
% sv_f001 = sv_f001/100;
% figure
% hold on
% plot((sv_f1), '.', 'MarkerSize', 20)
% plot((sv_f05), '.', 'MarkerSize', 20)
% plot((sv_f025), '.', 'MarkerSize', 20)
% plot((sv_f01), '.', 'MarkerSize', 20)
% plot((sv_f001), '.', 'MarkerSize', 20)
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Singular value")
% legend('.1', '.05', '0.25', '.01', '.001')
% title("Singular values of X for different frequency differences")
% hold off

%% Estimation of directions, Q2 - Noiseless check
% [M, N, Delta, theta, f, SNR] = deal(5, 20, 0.5, [-20; 30], [0.1; 0.3], 160);
% [X, A, S] = genData(M, N, Delta, theta, f, SNR);
% theta = esprit(X, 2)

%% Estimation of frequencies, Q3 - Noiseless check
% [M, N, Delta, theta, f, SNR] = deal(5, 20, 0.5, [-20; 30], [0.1; 0.3], 160);
% [X, A, S] = genData(M, N, Delta, theta, f, SNR);
% freq = espritfreq(X, 2)

%% Joint estimation of dirextions and frequencies, Q3 - Noiseless check
% [M, N, Delta, theta, f, SNR] = deal(5,20, 0.5, [-20; 30], [0.1; 0.3], 160);
% [X, A, S] = genData(M, N, Delta, theta, f, SNR);
% [theta, freq] = joint(X, 2, 5)

%% Comparison, Q1 - Plot of the estimation performance of the three algorithms
% avg_AoA_list = zeros(2, 11);
% avg_freq_list = zeros(1, 11);
% SNR_list = linspace(0,20,11)
% 
% % AoA_list = [];
% % freq_list = [];
% 
% AoA_list = zeros(2,11,1000);
% freq_list = zeros(2,11,1000);
%     
% % Take 1000 test runs 
% for j = 1:1000
%     % Run for each SNR value
%     for i = 1:11
%         SNR_instance = SNR_list(i);
%         [M, N, Delta, theta, f, SNR] = deal(3, 20, 0.5, [-20; 30], [0.1; 0.12], SNR_instance);
%         [X, A, S] = genData(M, N, Delta, theta, f, SNR);
%         AoA_list(:,i,j) = sort(esprit(X,2));
%         freq_list(:,i,j) = sort(espritfreq(X, 2));
%     end
% end
% 
% avg_AoA_list = mean(AoA_list,3);
% avg_freq_list = mean(freq_list,3);
% 
% std_AoA_list = std(squeeze(AoA_list(1,:,:)),0,2)
% std_freq_list = std(squeeze(freq_list(1,:,:)),0,2)
% 
% %plot AoA
% figure
% title("Mean values angle estimations")
% hold on
% plot(SNR_list, avg_AoA_list)
% 
% %plot freq
% figure
% title("Mean values frequency estimations")
% hold on
% plot(SNR_list, avg_freq_list)
% 
% %plot std AoA
% figure
% title("Std values angle estimations")
% hold on
% plot(SNR_list, std_AoA_list)
% 
% %plot std freq
% figure
% title("Std values frequency estimations")
% hold on
% plot(SNR_list, std_freq_list)

%% Comparison, Q2 - Compute zero-forcing beamformers, noiseless check
% % Generate data using the genData function
% [M, N, Delta, theta, f, SNR] = deal(5, 20, 0.5, [-20; 30], [0.1; 0.12], 160);
% [X, A, S] = genData(M, N, Delta, theta, f, SNR);
% % Using data estimate DoAs
% DoA = sort(esprit(X,2))
% % Use estimates DoAs to construct A
% phi = 2*pi*Delta*sin(deg2rad(DoA));
% A_hat = exp(1i*(0:M-1)'.*phi');
% W = pinv(A_hat);
% % Apply beamformer to obtain an estimate of S
% S_hat = W*X;
% round(S_hat, 5) == round(S, 5)
% % Use same data to estimate source frequencies
% freq = sort(espritfreq(X,2));
% S_hat_freq =  exp(1i*2*pi*freq*(0:N-1));
% A_hat = X*S_hat_freq'*inv(S_hat_freq*S_hat_freq');
% W_herm= pinv(A_hat);
% round(A_hat, 5) == round(A, 5)

%% Comparison, Q3 - Plot spatial response of beamformers
% [M, N, Delta, theta, f, SNR] = deal(3, 20, 0.5, [-20; 30], [0.1; 0.12], 10);
% [X, A, S] = genData(M, N, Delta, theta, f, SNR);
% % Using data estimate DoAs
% DoA = sort(esprit(X,2))
% % Use estimates DoAs to construct A
% phi = 2*pi*Delta*sin(deg2rad(DoA));
% A_hat = exp(1i*(0:M-1)'.*phi');
% W = pinv(A_hat);
% % Apply beamformer to obtain an estimate of S
% S_hat = W*X;
% round(S_hat, 5) == round(S, 5)
% % Use same data to estimate source frequencies
% freq = sort(espritfreq(X,2));
% S_hat_freq =  exp(1i*2*pi*freq*(0:N-1));
% A_hat = X*S_hat_freq'*inv(S_hat_freq*S_hat_freq');
% W_herm= pinv(A_hat);
% round(A_hat, 5) == round(A, 5)
% 
% % Plot spatial response for beamformer with A is known
% i=1
% for theta=-90:90
%     exponent = 2*pi*Delta*sin(deg2rad(theta));
%     a = exp(1i*(0:M-1)'.*exponent);
%     y(i, :) = abs(W*a);
%     i = i+1;
% end
% 
% figure
% title("DoA esprit")
% hold on
% plot(-90:90, y(:, 1))
% plot(-90:90, y(:, 2))
% xlabel("Theta(degree)")
% ylabel("|y(theta)|")
% hold off
% 
% % Plot spatial response for beamformer with S is known
% j=1
% for theta=-90:90
%     exponent = 2*pi*Delta*sin(deg2rad(theta));
%     a = exp(1i*(0:M-1)'.*exponent);
%     y_freq(j, :) = abs(W_herm*a);
%     j = j+1;
% end
% 
% figure
% title("Frequency Esprit")
% hold on
% plot(-90:90, y_freq(:, 1))
% plot(-90:90, y_freq(:, 2))
% xlabel("Theta(degree)")
% ylabel("|y(theta)|")
% hold off

%% Zero-forcing and Wiener equalizer, Q1 - Plot estimated symbols P=4
% N = 500;
% s = ones([N 1]);
% h = [1; -1; 1; -1];
% qpsk_symbols = [-1-1i, 1-1i, -1+1i, 1+1i]./sqrt(2);
% 
% for i=1:N
%     s(i) = qpsk_symbols(randi([1 4]));
% end
% 
% P = 4;
% sigma = 0.5;
% X = gendata_conv(s,P,N,sigma);
% 
% [S, s_hat, s_raw] = zeroForcingReceiver(h, X, N);
% [S_wf, s_hat_wf, s_raw_wf] = WienerReceiver(X, h,sigma, N);
% 
% % Plot estimated symbols of zero-forcing receiver
% figure;
% hold on;
% plot(s_raw_wf, '.', 'color', 'blue', 'MarkerSize', 5);
% plot(qpsk_symbols, 'x', 'color', 'red', 'MarkerSize', 15, 'LineWidth', 3);
% title("Raw estimates of Wiener receiver (N=500, P=4)")
% xlabel("Re(s)")
% ylabel("Im(s)")
% xlim([-1.5,  1.5]);
% ylim([-1.5,  1.5]);
% hold off;
% 
% % Plot estimated symbols of Wiener equalizer
% figure
% hold on;
% plot(s_raw, '.', 'color', 'blue', 'MarkerSize', 5);
% plot(qpsk_symbols, 'x', 'color', 'red', 'MarkerSize', 15, 'LineWidth', 3);
% title("Raw estimates of Zero-Forcing receiver (N=500, P=4)")
% xlabel("Re(s)")
% ylabel("Im(s)")
% xlim([-1.5,  1.5]);
% ylim([-1.5,  1.5]);
% hold off;

%% Zero-forcing and Wiener equalizer, Q2 - Plot estimated symbols P=8
% N = 500;
% s = ones([N 1]);
% h = [1; 1; -1; -1; 1; 1; -1; -1];
% qpsk_symbols = [-1-1i, 1-1i, -1+1i, 1+1i]./sqrt(2);
% 
% for i=1:N
%     s(i) = qpsk_symbols(randi([1 4]));
% end
% 
% P = 8;
% sigma = 0.5;
% X = gendata_conv(s,P,N,sigma);
% 
% [S, s_hat, s_raw] = zeroForcingReceiver(h, X, N);
% [S_wf, s_hat_wf, s_raw_wf] = WienerReceiver(X, h,sigma, N);
% 
% % Plot estimated symbols of zero-forcing receiver
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
% % Plot estimated symbols of Wiener equalizer
% figure
% hold on;
% plot(s_raw, '.', 'color', 'blue', 'MarkerSize', 5);
% plot(qpsk_symbols, 'x', 'color', 'red', 'MarkerSize', 15, 'LineWidth', 3);
% title("Raw estimates of Zero-Forcing receiver (N=500, P=8)")
% xlabel("Re(s)")
% ylabel("Im(s)")
% xlim([-1.5,  1.5]);
% ylim([-1.5,  1.5]);
% hold off;








% %% Q3, implement ESPRIT algorithm
% % Generate data using the genData function
[M, N, Delta, theta, f, SNR] = deal(3, 20, 0.5, [-20; 30], [0.1, 0.12]', 10);
[X, A, S] = genData(M, N, Delta, theta, f, SNR);
% % ESPRIT_freq(X, 2)
% %% Zero-forcing beamformer with A is known 
% 
% % % Obtain estimates of source directions, sort such that signal estimate
% % % will have the same order
% 
%BEAMFORMING ESPRIT DoA
DoA = sort(esprit(X,2));
% Construct zero-forcing beamformer
phi = 2*pi*Delta*sin(deg2rad(DoA));
A_hat = exp(1i*(0:M-1)'.*phi');
W = pinv(A_hat);
% Apply beamformer to obtain an estimate of S
S_hat = W*X;

%% Zero-forcing beamformer with S is known

% Obtain estimates of source directions, sort such that signal estimate
% will have the same order
% Construct zero-forcing beamformer
% Apply beamformer to obtain an estimate of S
A_hat = X*pinv(S_hat);      

i=1
for theta=-90:90
    exponent = 2*pi*Delta*sin(deg2rad(theta));
    a = exp(1i*(0:M-1)'.*exponent);
    y(i, :) = abs(W*a);
    i = i+1;
end

figure
title("DoA esprit")
hold on
plot(-90:90, y(:, 1))
plot(-90:90, y(:, 2))
xlabel("Theta(degree)")
ylabel("|y(theta)|")
hold off


freq = sort(espritfreq(X,2));
S_hat_freq =  exp(1i*2*pi*freq*(0:N-1));
A_hat = X*S_hat_freq'*inv(S_hat_freq*S_hat_freq');
W_herm= pinv(A_hat);
% W_herm = S_hat_freq*X'*inv(X*X')

j=1
for theta=-90:90
    exponent = 2*pi*Delta*sin(deg2rad(theta));
    a = exp(1i*(0:M-1)'.*exponent);
    y_freq(j, :) = abs(W_herm*a);
    j = j+1;
end

figure
title("Frequency Esprit")
hold on
plot(-90:90, y_freq(:, 1))
plot(-90:90, y_freq(:, 2))
xlabel("Theta(degree)")
ylabel("|y(theta)|")
hold off

%% Make a plot of the estimation performance of the three algorithms (mean values and
%% standard deviations of the estimated angles and frequencies)

% 
% avg_AoA_list = zeros(2, 11);
% avg_freq_list = zeros(1, 11);
% SNR_list = linspace(0,20,11);
% 
% 
% AoA_list = [];
% freq_list = [];
% 
%  AoA_list = zeros(2,11,1000);
%  freq_list = zeros(2,11,1000);
% 
% 
% for j = 1:1000
% 
%     for i = 1:11
%         SNR_instance = SNR_list(i);
%         [M, N, Delta, theta, f, SNR] = deal(3, 20, 0.5, [-20; 30], [0.2, 0.3]', SNR_instance);
%         [X, A, S] = genData(M, N, Delta, theta, f, SNR);
%         % AoA = sort(esprit(X,2));
%         % freq = espritfreq(X, 2);
%         AoA_list(:,i,j) = sort(esprit(X,2));
%         freq_list(:,i,j) = sort(espritfreq(X, 2));
% 
% 
%         % AoA_list = [AoA_list, AoA];
%         % freq_list = [freq_list, freq];
%     end
% 
% 
%     % avg_AoA_list = avg_AoA_list + AoA_list;
%     % avg_freq_list = avg_freq_list + freq_list;
% end
% 
% avg_AoA_list = mean(AoA_list,3);
% avg_freq_list = mean(freq_list,3);
% 
% 
% std_AoA_list = std(squeeze(AoA_list(1,:,:)),0,2)
% std_freq_list = std(squeeze(freq_list(1,:,:)),0,2)

% std_freq_list = std(freq_list);
% avg_AoA_list = avg_AoA_list ./ 1000;
% avg_freq_list = avg_freq_list ./ 1000;

%plot AoA
% figure
% title("Mean values angle estimations")
% hold on
% plot(SNR_list,avg_AoA_list)

%plot freq
% figure
% title("Mean values frequency estimations")
% hold on
% plot(SNR_list,avg_freq_list)

%plot std AoA
% figure
% title("Std values angle estimations")
% hold on
% plot(SNR_list,std_AoA_list)

%plot std freq
% figure
% title("Std values frequency estimations")
% hold on
% plot(SNR_list,std_freq_list)



% N = 500;
% s = ones([N 1]);
% 
% qpsk_symbols = [-1-1i, 1-1i, -1+1i, 1+1i]./sqrt(2);
% 
% for i=1:N
%     s(i) = qpsk_symbols(randi([1 4]));
% end
% 
% P = 4;
% sigma = 0.0;
% X = gendata_conv(s,P,N,sigma);
% % x = gendata_conv_2(s,P,N,sigma);
% rank(X)
% 
% h = [1; -1; 1; -1];
% nul = zeros([4 1]);
% H = [h, nul; nul, h]
% S_hat = pinv(H)*X;
% 
% 
% %% Wiener receiver
% 
%

% S_wiener = WienerReceiver(X,h);

%
%
%
%
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





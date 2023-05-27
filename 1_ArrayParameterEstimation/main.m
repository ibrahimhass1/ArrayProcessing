%% Q3, implement ESPRIT algorithm
% Generate data using the genData function
[M, N, Delta, theta, f, SNR] = deal(5, 20, 0.5, [-20; 30], [0.1; 0.3], 20);
[X, A, S] = genDataKR(M, N, Delta, theta, f, SNR);
% [X, A, S] = genDataKR(M, N, Delta, theta, f, SNR);
% size(X)
% % Call the esprit function
% AoA = esprit(X, 2)
% freq = espritfreq(X, 2)

% Apply temporal smoothing, with a factor m
X = eye(4)
m = 1;
Xs = temporal_smoothing(X,m)

%% Q2, effect of number of samples on the singular values
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
% plot(log10(sv_N100), '.', 'MarkerSize', 20)
% plot(log10(sv_N60), '.', 'MarkerSize', 20)
% plot(log10(sv_N40), '.', 'MarkerSize', 20)
% 
% plot(log10(sv_N20), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Log10 of singular value")
% legend('N=160', 'N=80', 'N=40', 'N=20')
% title("Singular values for different number of samples")
% hold off
% 
% 10^mean(log10(sv_N100)-log10(sv_N20))
% 10^mean(log10(sv_N60)-log10(sv_N20))
% 10^mean(log10(sv_N40)-log10(sv_N20))

%% Q2, effect of number of antennas on the singular values
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
% plot(log(sv_A15), '.', 'MarkerSize', 20)
% plot(log(sv_A10), '.', 'MarkerSize', 20)
% plot(log(sv_A5), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:15) 
% ylabel("Log10 of singular value")
% legend('M=15', 'M=10', 'M=5')
% title("Singular values for different number of antennas")
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
% plot(log10(sv_A50), '.', 'MarkerSize', 20)
% plot(log10(sv_A10), '.', 'MarkerSize', 20)
% plot(log10(sv_A1), '.', 'MarkerSize', 20)
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Log10 of singular value")
% legend('50 deg', '10 deg', '1 deg')
% title("Singular values for different angles between sources")
% hold off
% 

%% Q2, effect of the source frequency difference on the singular values
% sv_f2 = 0;
% sv_f1 = 0;
% sv_f05 = 0;
% sv_f01 = 0;
% sv_f001 = 0;
% 
% for i=1:100
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.1;.3],20);
%     sv_f2 =sv_f2 +svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.2;.3],20);
%     sv_f1=sv_f1+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.25;.3],20);
%     sv_f05=sv_f05+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.29;.3],20);
%     sv_f01=sv_f01+svd(X);
%     [X,A,S] = genData(5,20,.5,[-20; 30],[.29;.3],20);
%     sv_f001=sv_f001+svd(X);
% end
% sv_f2 = sv_f2/100
% sv_f1 = sv_f1/100
% sv_f05 = sv_f05/100
% sv_f01 = sv_f01/100
% 
% 
% figure
% hold on
% plot(log10(sv_f2), '.', 'MarkerSize', 20)
% plot(log10(sv_f1), '.', 'MarkerSize', 20)
% plot(log10(sv_f05), '.', 'MarkerSize', 20)
% plot(log10(sv_f01), '.', 'MarkerSize', 20)
% 
% 
% xlabel("Singular value index")
% xticks(1:5) 
% ylabel("Log10 of singular value")
% legend('.2', '.1', '.05', '.01')
% title("Singular values for different frequency differences")
% hold off





load('snr_values.mat')

nanIndices = isnan(SNR_MF_09);
snr_mean_wf_90= SNR_MF_09(~nanIndices);

nanIndices = isnan(SNR_MF_05);
snr_mean_wf_50= SNR_MF_05(~nanIndices);

nanIndices = isnan(SNR_MF_01);
snr_mean_wf_10= SNR_MF_01(~nanIndices);

snr_mean_wf_90 = mean(snr_mean_wf_90)
snr_mean_wf_50 = mean(snr_mean_wf_50)
snr_mean_wf_10 = mean(snr_mean_wf_10)
wf = [snr_mean_wf_10 snr_mean_wf_50 snr_mean_wf_90]
wf = 10*log10(wf)
wf_mse = [MSE_WF_01 MSE_WF_05 MSE_WF_09]

nanIndices = isnan(SNR_DaS_09);
snr_mean_DaS_90= SNR_DaS_09(~nanIndices);

nanIndices = isnan(SNR_DaS_05);
snr_mean_DaS_50= SNR_DaS_05(~nanIndices);

nanIndices = isnan(SNR_DaS_01);
snr_mean_DaS_10= SNR_DaS_01(~nanIndices);

snr_mean_DaS_90 = mean(snr_mean_DaS_90)
snr_mean_DaS_50 = mean(snr_mean_DaS_50)
snr_mean_DaS_10 = mean(snr_mean_DaS_10)

das = [snr_mean_DaS_10 snr_mean_DaS_50 snr_mean_DaS_90]
das = 10*log10(das)

das_mse = [MSE_DaS_01 MSE_DaS_05 MSE_DaS_09]

nanIndices = isnan(SNR_MVDR_09);
snr_mean_MVDR_90= SNR_MVDR_09(~nanIndices);

nanIndices = isnan(SNR_MVDR_05);
snr_mean_MVDR_50= SNR_MVDR_05(~nanIndices);

nanIndices = isnan(SNR_MVDR_01);
snr_mean_MVDR_10= SNR_MVDR_01(~nanIndices);

snr_mean_MVDR_90 = mean(snr_mean_MVDR_90)
snr_mean_MVDR_50 = mean(snr_mean_MVDR_50)
snr_mean_MVDR_10 = mean(snr_mean_MVDR_10)
mvdr = [snr_mean_MVDR_10 snr_mean_MVDR_50 snr_mean_MVDR_90]
mvdr = 10*log10(mvdr)
mvdr_mse = [MSE_MVDR_01 MSE_MVDR_05 MSE_MVDR_09]

% 

figure('Position', [100, 100, 800, 600]);
% Bar chart 1
subplot(2, 3, 1);
x = [1 2 3]
bar(x, das);
% Define the labels for each bar
labels = {'0.1', '0.5', '0.9'};

% Set the x-axis tick labels
xticklabels(labels);
ylabel("SNR[dB]")

ylim([0,  4]);
title('Delay-and-Sum');
hold on;

% Bar chart 3
subplot(2, 3, 2);
bar(x,mvdr);
% Set the x-axis tick labels
xticklabels(labels);
ylabel("SNR[dB]")

ylim([0,  4]);
title('MVDR');
hold on;

% Bar chart 5
subplot(2, 3, 3);
bar(x, wf);
ylabel("SNR[dB]")
% Set the x-axis tick labels
xticklabels(labels);
ylim([0,  4]);
title('Multi-Channel Wiener');
hold on;


subplot(2, 3, 4);
bar(x,das_mse);
ylabel("MSE")
% Set the x-axis tick labels
xticklabels(labels);
xlabel("alpha")
ylim([0,  0.004]);
title('Delay-and-Sum');
hold on;

% Bar chart 3
subplot(2, 3, 5);
bar(x,mvdr_mse);
ylabel("MSE")
% Set the x-axis tick labels
xticklabels(labels);
xlabel("alpha")
ylim([0,  0.004]);
title('MVDR');
hold on;

% Bar chart 5
subplot(2, 3, 6);
bar(x, wf_mse);
ylabel("MSE")
% Set the x-axis tick labels
xticklabels(labels);
xlabel("alpha")
ylim([0,  0.005]);
title('Multi-Channel Wiener');
hold on;

% % Bar chart 6
% subplot(2, 2, 4);
% bar(values6);
% ylabel("Gain")
% xlabel("Number of users")
% ylim([.85,  1.75]);
% title('Proportional Fair SRB-basis');
% hold on;
% % Set color for bars lower than 1 to red
% bar(find(values6 < 1), values6(values6 < 1), 'r');
% 
% % Adjusting the figure layout
sgtitle('SNR and MSE for the different beamformers, white noise case');

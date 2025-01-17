
load('snr_values.mat')
nanIndices = isnan(SNR_MF_09_NW);
snr_mean_wf_90_NW= SNR_MF_09_NW(~nanIndices);

nanIndices = isnan(SNR_MF_05_NW);
snr_mean_wf_50_NW= SNR_MF_05_NW(~nanIndices);

nanIndices = isnan(SNR_MF_01_NW);
snr_mean_wf_10_NW= SNR_MF_01_NW(~nanIndices);

snr_mean_wf_90_NW = mean(snr_mean_wf_90_NW)
snr_mean_wf_50_NW = mean(snr_mean_wf_50_NW)
snr_mean_wf_10_NW = mean(snr_mean_wf_10_NW)
wf_nw = [snr_mean_wf_10_NW snr_mean_wf_50_NW snr_mean_wf_90_NW]
wf_nw = 10*log10(wf_nw)


nanIndices = isnan(SNR_DaS_09_NW);
snr_mean_DaS_90_NW= SNR_DaS_09_NW(~nanIndices);

nanIndices = isnan(SNR_DaS_05_NW);
snr_mean_DaS_50_NW= SNR_DaS_05_NW(~nanIndices);

nanIndices = isnan(SNR_DaS_01_NW);
snr_mean_DaS_10_NW= SNR_DaS_01_NW(~nanIndices);

snr_mean_DaS_90_NW = mean(snr_mean_DaS_90_NW)
snr_mean_DaS_50_NW = mean(snr_mean_DaS_50_NW)
snr_mean_DaS_10_NW = mean(snr_mean_DaS_10_NW)

das_nw = [snr_mean_DaS_10_NW snr_mean_DaS_50_NW snr_mean_DaS_90_NW]
das_nw = 10*log10(das_nw)

das_mse = [MSE_DaS_01 MSE_DaS_05 MSE_DaS_09]

nanIndices = isnan(SNR_MVDR_09_NW);
snr_mean_MVDR_90_NW= SNR_MVDR_09_NW(~nanIndices);

nanIndices = isnan(SNR_MVDR_05_NW);
snr_mean_MVDR_50_NW= SNR_MVDR_05_NW(~nanIndices);

nanIndices = isnan(SNR_MVDR_01_NW);
snr_mean_MVDR_10_NW= SNR_MVDR_01_NW(~nanIndices);

snr_mean_MVDR_90_NW = mean(snr_mean_MVDR_90_NW)
snr_mean_MVDR_50_NW = mean(snr_mean_MVDR_50_NW)
snr_mean_MVDR_10_NW = mean(snr_mean_MVDR_10_NW)
mvdr_nw = [snr_mean_MVDR_10_NW snr_mean_MVDR_50_NW snr_mean_MVDR_90_NW]
mvdr_nw = 10*log10(mvdr_nw)
mvdr_mse = [MSE_MVDR_01 MSE_MVDR_05 MSE_MVDR_09]

% 

figure('Position', [100, 100, 800, 300]);
% Bar chart 1
subplot(1, 3, 1);
x = [1 2 3]
bar(x,das);
hold on
bar(x, das_nw);
% Define the labels for each bar
labels = {'0.1', '0.5', '0.9'};

% Set the x-axis tick labels
xticklabels(labels);
ylabel("SNR[dB]")
xlabel("alpha")

ylim([0,  4]);
title('Delay-and-Sum');
hold on;

% Bar chart 3
subplot(1, 3, 2);

bar(x,mvdr);
hold on
bar(x,mvdr_nw);
% Set the x-axis tick labels
xticklabels(labels);
ylabel("SNR[dB]")
xlabel("alpha")

ylim([0,  4]);
title('MVDR');
hold on;

% Bar chart 5
subplot(1, 3, 3);
bar(x,wf);
hold on
bar(x, wf_nw);
legend('White', 'Non-white');
ylabel("SNR[dB]")
xlabel("alpha")
% Set the x-axis tick labels
xticklabels(labels);
ylim([0,  4]);
title('Multi-Channel Wiener');
hold on;

sgtitle('SNR for the different beamformers');

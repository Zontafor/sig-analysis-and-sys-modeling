%% hw6.m
%  White-noise driven LTI system
%  Correlation, spectra, TF and IRF estimation

clear; close all; clc;

%% Define paths

data_path = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw6/data';
fig_path = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw6/figs';

%% Save fig function

function save_fig(fh, fname)
    set(fh,'Units','inches');
    pos = get(fh,'Position');
    set(fh,'PaperPositionMode','auto');
    set(fh,'PaperUnits','inches');
    set(fh,'PaperPosition',[0 0 pos(3) pos(4)]);
    set(fh,'PaperSize',[pos(3) pos(4)]);
    print(fh, fname, '-dpdf', '-r300');
end

%% Parameters and true IRF

M       = 50;
maxLag  = 2*M;
N_short = 4000;
N_long  = 16000;

m      = 0:M-1;
h_true = 2*(exp(-m/8) - exp(-m/4)).';   % true IRF, column

Nfft = 2^nextpow2(N_long + M - 1);

corr_est = @(x,y,L) deal(xcorr(x,y,L,'unbiased'), -L:L);

%% Data generation and correlation estimates

rng(1);
w_short = randn(N_short,1);
y_short = filter(h_true,1,w_short);

rng(2);
w_long  = randn(N_long,1);
y_long  = filter(h_true,1,w_long);

[Rww_s,lags_s] = corr_est(w_short, w_short, maxLag);
[Ryy_s,~]      = corr_est(y_short, y_short, maxLag);
[Rwy_s,~]      = corr_est(w_short, y_short, maxLag);

[Rww_l,lags_l] = corr_est(w_long, w_long, maxLag);
[Ryy_l,~]      = corr_est(y_long, y_long, maxLag);
[Rwy_l,~]      = corr_est(w_long, y_long, maxLag);

%% Direct time-domain IRF estimate from (w_long, y_long)

sigma2_hat = var(w_long,1);
h_td_est_l = zeros(M,1);

for k = 0:M-1
    num = w_long(1:N_long-k).' * y_long(1+k:N_long);
    h_td_est_l(k+1) = num / ((N_long-k) * sigma2_hat);
end

%% Spectra, cross-spectra, and TF estimates

W_s = fft(w_short, Nfft);
Y_s = fft(y_short, Nfft);
Sw_s  = (1/N_short) * (W_s .* conj(W_s));
Swy_s = (1/N_short) * (conj(W_s) .* Y_s);
Hhat_s = Swy_s ./ (Sw_s + eps);

W_l = fft(w_long, Nfft);
Y_l = fft(y_long, Nfft);
Sw_l  = (1/N_long) * (W_l .* conj(W_l));
Swy_l = (1/N_long) * (conj(W_l) .* Y_l);
Hhat_l = Swy_l ./ (Sw_l + eps);

H_true = fft(h_true, Nfft);

%% IRF estimate via inverse DFT (long record)

h_fd_est_l = real(ifft(Hhat_l));
h_fd_est_l = h_fd_est_l(1:M);

%% Model prediction and NMSE on short record

p_td_short = filter(h_td_est_l, 1, w_short);
p_fd_short = filter(h_fd_est_l, 1, w_short);

num_td = mean((p_td_short - y_short).^2);
num_fd = mean((p_fd_short - y_short).^2);
den_y  = mean(y_short.^2);

NMSE_td = num_td / den_y;
NMSE_fd = num_fd / den_y;

fprintf('NMSE (time-domain IRF, short data)      = %.6g\n', NMSE_td);
fprintf('NMSE (frequency-domain IRF, short data) = %.6g\n', NMSE_fd);

%% Least-squares IRF estimate (time domain, long record)

L    = N_long - M + 1;
Wmat = zeros(L, M);

for n = M:N_long
    Wmat(n-M+1,:) = w_long(n:-1:n-M+1).';
end
y_vec = y_long(M:N_long);

h_ls = (Wmat' * Wmat) \ (Wmat' * y_vec);

rel_err = @(h_est) norm(h_est - h_true) / norm(h_true);

err_td = rel_err(h_td_est_l);
err_fd = rel_err(h_fd_est_l);
err_ls = rel_err(h_ls);

fprintf('Relative IRF error (time-domain, direct) = %.6g\n', err_td);
fprintf('Relative IRF error (frequency-domain TF) = %.6g\n', err_fd);
fprintf('Relative IRF error (LS time-domain)      = %.6g\n', err_ls);

%% Plotting correlations for N = 4000 vs 16000

fig_corr = figure;
subplot(3,1,1);
plot(lags_s, Rww_s, 'b-', lags_l, Rww_l, 'r-');
xlabel('Lag k'); ylabel('R_{ww}(k)');
legend('N = 4000','N = 16000','Location','best');
title('Input autocorrelation');

subplot(3,1,2);
plot(lags_s, Ryy_s, 'b-', lags_l, Ryy_l, 'r-');
xlabel('Lag k'); ylabel('R_{yy}(k)');
legend('N = 4000','N = 16000','Location','best');
title('Output autocorrelation');

subplot(3,1,3);
plot(lags_s, Rwy_s, 'b-', lags_l, Rwy_l, 'r-');
xlabel('Lag k'); ylabel('R_{wy}(k)');
legend('N = 4000','N = 16000','Location','best');
title('Input–output cross-correlation');

set(fig_corr, 'Units','inches','Position',[1 1 6 6]);
save_fig(fig_corr, fullfile(fig_path,'hw6_fig1_corr_ww_yy_wy'));

%% Plotting magnitude and phase of TF estimates vs true TF

k_max = floor(Nfft/2);                % plot up to Nyquist
omega = (0:k_max-1)/k_max * pi;       % normalized rad/sample

H_true_half = H_true(1:k_max);
Hhat_s_half = Hhat_s(1:k_max);
Hhat_l_half = Hhat_l(1:k_max);

fig_tf = figure;
subplot(2,1,1);
plot(omega, 20*log10(abs(H_true_half)), 'k-', ...
     omega, 20*log10(abs(Hhat_s_half)), 'b--', ...
     omega, 20*log10(abs(Hhat_l_half)), 'r-.');
xlabel('\omega (rad/sample)'); ylabel('|H(e^{j\omega})| (dB)');
legend('True','N = 4000','N = 16000','Location','best');
title('Magnitude response');

subplot(2,1,2);
plot(omega, unwrap(angle(H_true_half))*180/pi, 'k-', ...
     omega, unwrap(angle(Hhat_s_half))*180/pi, 'b--', ...
     omega, unwrap(angle(Hhat_l_half))*180/pi, 'r-.');
xlabel('\omega (rad/sample)'); ylabel('Phase (deg)');
legend('True','N = 4000','N = 16000','Location','best');
title('Phase response');

set(fig_tf, 'Units','inches','Position',[1 1 6 6]);
save_fig(fig_tf, fullfile(fig_path,'hw6_fig2_tf_mag_phase'));

%% Plotting IRF estimates vs true IRF

n_irf = 0:M-1;

fig_irf = figure;
stem(n_irf, h_true, 'k','filled'); hold on;
stem(n_irf+0.05, h_td_est_l, 'b');    % small shifts for visibility
stem(n_irf+0.10, h_fd_est_l, 'r');
stem(n_irf+0.15, h_ls,        'g');
hold off;
xlabel('m'); ylabel('h[m]');
legend('True','Time-domain (corr)','Freq-domain (IDFT)','LS', ...
       'Location','northeast');
title('IRF estimates');

set(fig_irf, 'Units','inches','Position',[1 1 6 4]);
save_fig(fig_irf, fullfile(fig_path,'hw6_fig3_irf_estimates'));

%% Plotting output prediction on short record

n_view = 1:200;

fig_pred = figure;

% Colorblind palette
cb_blue   = [0 0.4470 0.7410];
cb_orange = [0.8500 0.3250 0.0980];
cb_green  = [0.4660 0.6740 0.1880];

plot(n_view, y_short(n_view), 'Color', 'k', 'LineWidth', 1.5); hold on;
plot(n_view, p_td_short(n_view), '--', 'Color', cb_blue, 'LineWidth', 1.5);
plot(n_view, p_fd_short(n_view), '-.', 'Color', cb_orange, 'LineWidth', 1.5);

xlabel('n'); ylabel('Amplitude');
legend('True y[n]', ...
       'Pred (time-domain IRF)', ...
       'Pred (freq-domain IRF)', ...
       'Location','best');
title('Output prediction on short record (N = 4000)');
hold off;

% % Normal palette
% plot(n_view, y_short(n_view), 'k-', ...
%      n_view, p_td_short(n_view), 'b--', ...
%      n_view, p_fd_short(n_view), 'r-.');
% xlabel('n'); ylabel('Amplitude');
% legend('True y[n]','Pred (time-domain IRF)','Pred (freq-domain IRF)', ...
%        'Location','best');
% title('Output prediction on short record (N = 4000)');

set(fig_pred, 'Units','inches','Position',[1 1 6 4]);
save_fig(fig_pred, fullfile(fig_path,'hw6_fig4_output_prediction'));

%% Optional saving for numerical outputs
%  .mat structure

results = struct;
results.NMSE_time_domain    = NMSE_td;
results.NMSE_freq_domain    = NMSE_fd;
results.relerr_time_domain  = err_td;
results.relerr_freq_domain  = err_fd;
results.relerr_LS           = err_ls;

save(fullfile(data_path, 'hw6_results.mat'), 'results');

% .txt or .csv
fid = fopen(fullfile(data_path, 'hw6_results_summary.txt'), 'w');

fprintf(fid, 'HW6 Numerical Results Summary\n');
fprintf(fid, '-----------------------------\n');
fprintf(fid, 'NMSE (time-domain IRF):      %.10f\n', NMSE_td);
fprintf(fid, 'NMSE (freq-domain IRF):      %.10f\n', NMSE_fd);
fprintf(fid, 'Relative error (TD IRF):     %.10f\n', err_td);
fprintf(fid, 'Relative error (FD IRF):     %.10f\n', err_fd);
fprintf(fid, 'Relative error (LS IRF):     %.10e\n', err_ls);

fclose(fid);

% CSV for Python or R
T = table(NMSE_td, NMSE_fd, err_td, err_fd, err_ls, ...
          'VariableNames', {'NMSE_TD','NMSE_FD','Err_TD','Err_FD','Err_LS'});

writetable(T, fullfile(data_path,'hw6_results.csv'));
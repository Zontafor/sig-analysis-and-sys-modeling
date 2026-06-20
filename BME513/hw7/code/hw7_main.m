%% hw7_main.m

clear; clc; close all;

% define paths and load data
data_dir = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw7/code/data';
fig_dir  = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw7/code/figs';

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

data_file = fullfile(data_dir, 'subj_data.mat');
S         = load(data_file);        % expects S.data, S.deltaT
data      = S.data;
Ts        = S.deltaT(1);
fs        = 1/Ts;

bp    = data.bp(:);
fv    = data.fv(:);
hr    = data.hr(:);
etco2 = data.etco2(:);

N      = numel(bp);
t      = (0:N-1).' * Ts;
maxlag = 100;

% Q1: Histograms + Gaussian fits + bp–fv scatter
vars  = {bp, fv, hr, etco2};
names = {'bp', 'fv', 'hr', 'etco2'};
nbins = 20;

fig_q1_hist = figure;
for k = 1:4
    x = vars{k};
    subplot(2,2,k);
    
    % histogram as PDF
    [counts, edges] = histcounts(x, nbins, 'Normalization', 'pdf');
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    bar(centers, counts, 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold on;
    
    % Gaussian PDF fit
    mu  = mean(x);
    sig = std(x);
    pdf_gauss = (1/(sig*sqrt(2*pi))) * exp(-0.5*((centers-mu)/sig).^2);
    plot(centers, pdf_gauss, 'LineWidth', 1.5);
    
    xlabel(sprintf('%s amplitude', names{k}));
    ylabel('PDF');
    title(sprintf('%s: hist + Gaussian fit', names{k}));
    grid on;
end

outfile = fullfile(fig_dir, 'HW7_Q1_histograms.png');
exportgraphics(fig_q1_hist, outfile, 'Resolution', 300);

fig_q1_scatter = figure;
plot(bp, fv, '.', 'MarkerSize', 6);
xlabel('bp (mmHg, demeaned)');
ylabel('fv (cm/s, demeaned)');
title('Contemporaneous bp vs fv');
grid on; axis tight;

outfile = fullfile(fig_dir, 'HW7_Q1_bp_fv_scatter.png');
exportgraphics(fig_q1_scatter, outfile, 'Resolution', 300);

% Q2: DFT magnitudes vs PSD from autocorrelation
Lac   = maxlag;
Nfft  = 2048;
faxis = (0:Nfft-1).' * (fs/Nfft);

fig_q2 = figure;
for k = 1:4
    x = vars{k};
    
    % raw DFT magnitude
    X    = fft(x, Nfft);
    magX = abs(X);
    
    % PSD from autocorrelation
    [rxx, ~] = xcorr(x, Lac, 'biased');
    Sxx      = real(fft(rxx, Nfft));
    
    subplot(2,2,k);
    plot(faxis, magX, 'LineWidth', 1); hold on;
    plot(faxis, Sxx,  'LineWidth', 1);
    xlim([0 fs/2]);
    xlabel('f (Hz)');
    ylabel('Amplitude / PSD (arb)');
    title(sprintf('%s: |X(f)| vs S_{xx}(f)', names{k}));
    legend('|X(f)|', 'S_{xx}(f) from r_{xx}', 'Location', 'best');
    grid on;
end

outfile = fullfile(fig_dir, 'HW7_Q2_DFT_vs_PSD.png');
exportgraphics(fig_q2, outfile, 'Resolution', 300);

% Q3: Transfer functions from cross-spectra
x_list = {bp,     bp,      etco2,   etco2};
y_list = {fv,     hr,      fv,      hr};
labels = {'bp->fv','bp->hr','etco2->fv','etco2->hr'};

Lxy       = maxlag;
Hhat_all  = cell(1,4);

fig_q3 = figure;
for k = 1:4
    x = x_list{k};
    y = y_list{k};
    
    % auto- and cross-correlations
    [rxx, ~] = xcorr(x, Lxy, 'biased');
    [ryx, ~] = xcorr(y, x, Lxy, 'biased');
    
    % spectra
    Sxx = fft(rxx, Nfft);
    Syx = fft(ryx, Nfft);
    
    % TF estimate
    Hhat        = Syx ./ (Sxx + eps);
    Hhat_all{k} = Hhat;
    
    magH = abs(Hhat);
    phH  = unwrap(angle(Hhat));
    
    subplot(4,2,2*k-1);
    plot(faxis, magH, 'LineWidth', 1);
    xlim([0 fs/2]);
    xlabel('f (Hz)');
    ylabel('|H(f)|');
    title(sprintf('%s: |H(f)|', labels{k}));
    grid on;
    
    subplot(4,2,2*k);
    plot(faxis, phH*180/pi, 'LineWidth', 1);
    xlim([0 fs/2]);
    xlabel('f (Hz)');
    ylabel('\angle H(f) (deg)');
    title(sprintf('%s: \angle H(f)', labels{k}));
    grid on;
end

outfile = fullfile(fig_dir, 'HW7_Q3_transfer_functions.png');
exportgraphics(fig_q3, outfile, 'Resolution', 300);

% Q4: IRFs from inverse TF and LS FIR
Lirf  = 40;
t_irf = (0:Lirf-1).' * Ts;

irf_ifft = cell(1,4);
irf_ls   = cell(1,4);

fig_q4 = figure;
for k = 1:4
    x    = x_list{k};
    y    = y_list{k};
    Hhat = Hhat_all{k};
    
    % IRF from inverse FFT
    h_circ      = real(ifft(Hhat, Nfft));
    h_ifft      = h_circ(1:Lirf);
    irf_ifft{k} = h_ifft;
    
    % LS FIR model
    n0   = Lirf;
    nmax = numel(x)-1;
    Mls  = nmax - n0 + 1;
    Xreg = zeros(Mls, Lirf);
    yreg = zeros(Mls, 1);
    for i = 1:Mls
        n = n0 + i - 1;
        Xreg(i,:) = flipud(x(n-Lirf+1:n));
        yreg(i)   = y(n);
    end
    h_ls      = Xreg \ yreg;
    irf_ls{k} = h_ls;
    
    subplot(2,2,k);
    plot(t_irf, h_ifft, 'LineWidth', 1.5); hold on;
    plot(t_irf, h_ls,   '--', 'LineWidth', 1.5);
    xlabel('time lag (s)');
    ylabel('h(t)');
    title(sprintf('IRF: %s', labels{k}));
    legend('IFFT(H)','LS-FIR','Location','best');
    grid on;
end

outfile = fullfile(fig_dir, 'HW7_Q4_IRFs_IFFT_vs_LS.png');
exportgraphics(fig_q4, outfile, 'Resolution', 300);

%  Q5: IRFs via Laguerre Expansion Technique (LET-1)
Llag  = 4;      % number of Laguerre functions
Mlag  = 50;     % memory length
alpha = 0.5;    % Laguerre parameter (a = sqrt(alpha))

h_let = cell(1,4);
Q_let = zeros(4,1);
t_let = (0:Mlag-1).' * Ts;

fig_q5 = figure;
for k = 1:4
    x = x_list{k};
    y = y_list{k};
    
    [hest, yest, Gest, Q, LC] = LET1_single_pair(x, y, Llag, Mlag, alpha);
    h_let{k} = hest;
    Q_let(k) = Q;
    
    subplot(2,2,k);
    plot(t_let, hest, 'LineWidth', 1.5);
    xlabel('time lag (s)');
    ylabel('h_{LET}(t)');
    title(sprintf('LET IRF: %s (Q=%.3f)', labels{k}, Q));
    grid on;
end

outfile = fullfile(fig_dir, 'HW7_Q5_IRFs_LET.png');
exportgraphics(fig_q5, outfile, 'Resolution', 300);

% optional: save workspace snapshot
save(fullfile(fig_dir, 'HW7_results_workspace.mat'));
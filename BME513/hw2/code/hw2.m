clear; clc;

%% BME 513 - Homework 2: Problem 3
% Author: Michelle L. Wu
% 
% Objective: Generate discrete-time signals x(n) and y(n), compute their DFTs,
% and analyze the resulting spectra.

%% Paths
fig_dir  = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw2/figs';
data_dir = '/Users/mlwu/Documents/Academia/USC/BME/513/hw/hw2/data_out';
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% Parameters
N   = 1024;              % original record length
n   = 1:N;               % problem specifies n = 1..1024
ZP  = 8;                 % zero-padding factor for finer freq sampling
Nf  = ZP*N;              % FFT length (interpolates spectrum)
epsdB = eps;             % to avoid log of zero


%% Signal definitions (exact)
% x(n) = exp(-n/400) * [cos(n*pi/64) - sin(n*pi/8)]
x = exp(-n/400) .* (cos(n*pi/64) - sin(n*pi/8));

% y(n) = cos(pi*n^2/8192) + sin(n*pi/32)
y = cos((n.^2)*pi/8192) + sin(n*pi/32);

%% DFTs (centered, amplitude-normalized to record length N)
X = fftshift(fft(x, Nf))/N;
Y = fftshift(fft(y, Nf))/N;

%% Frequency axis (cycles/sample) in [-0.5, 0.5)
f = (-Nf/2:Nf/2-1)/Nf;

%% Magnitude/phase
magX = abs(X);  phX = unwrap(angle(X));
magY = abs(Y);  phY = unwrap(angle(Y));

%% Expected resonances (cycles/sample)
fx_peaks = [1/128, 1/16];   % x(n) components
fy_peaks = [1/64];          % y(n) pure tone
fy_chirp_band = [0, 1/8];   % y(n) chirp sweep

%% Compute -3 dB bandwidths (positive peaks only; symmetric at -f)
BWrows = struct('Signal',{}, 'f0',{}, 'f_left_m3dB',{}, 'f_right_m3dB',{}, ...
                'BW_m3dB',{}, 'Peak_dB',{}, 'Peak_freq',{}); %#ok<*STRNU>

k = 0;
for f0 = fx_peaks
    [fL, fR, BW, pkdB, fpk] = minus3dB_bw(f, magX, f0);
    k = k+1; BWrows(k) = struct('Signal',"x(n)", 'f0',f0, ...
        'f_left_m3dB',fL, 'f_right_m3dB',fR, 'BW_m3dB',BW, ...
        'Peak_dB',pkdB, 'Peak_freq',fpk);
end
for f0 = fy_peaks
    [fL, fR, BW, pkdB, fpk] = minus3dB_bw(f, magY, f0);
    k = k+1; BWrows(k) = struct('Signal',"y(n)", 'f0',f0, ...
        'f_left_m3dB',fL, 'f_right_m3dB',fR, 'BW_m3dB',BW, ...
        'Peak_dB',pkdB, 'Peak_freq',fpk);
end

bw_tbl = struct2table(BWrows);   % already numeric; no cell2mat needed
fprintf('\nMeasured -3 dB Bandwidths (cycles/sample)\n');
disp(bw_tbl);
csv_path = fullfile(data_dir, 'measured_bandwidths.csv');
writetable(bw_tbl, csv_path);

fig1 = figure('Color','w'); tiledlayout(2,1);

nexttile
plot(f, 20*log10(magX+epsdB), 'LineWidth', 1.2); grid on; hold on
xline(fx_peaks(1),'--','f=1/128'); xline(-fx_peaks(1),'--');
xline(fx_peaks(2),'--','f=1/16');  xline(-fx_peaks(2),'--');
for f0 = fx_peaks
    [fL,fR] = minus3dB_bw(f, magX, f0);
    xline(fL,':'); xline(fR,':');
end
xlabel('Normalized frequency (cycles/sample)');
ylabel('|X| [dB]');
title('Magnitude Spectrum of x(n)');
ylim([-100 10]); hold off

nexttile
plot(f, phX, 'LineWidth', 1.2); grid on
xlabel('Normalized frequency (cycles/sample)');
ylabel('Phase ∠X [rad]');
title('Phase Spectrum of x(n)');
saveas(fig1, fullfile(fig_dir, 'x_spectrum.png'));

fig2 = figure('Color','w'); tiledlayout(2,1);

nexttile
plot(f, 20*log10(magY+epsdB), 'LineWidth', 1.2); grid on; hold on
xline(fy_peaks(1),'--','f=1/64'); xline(-fy_peaks(1),'--');
xline(fy_chirp_band(1),':','Chirp Start'); xline(fy_chirp_band(2),':','Chirp End');
[fL,fR] = minus3dB_bw(f, magY, fy_peaks(1));
xline(fL,':'); xline(fR,':');
xlabel('Normalized frequency (cycles/sample)');
ylabel('|Y| [dB]');
title('Magnitude Spectrum of y(n)');
ylim([-100 10]); hold off

nexttile
plot(f, phY, 'LineWidth', 1.2); grid on
xlabel('Normalized frequency (cycles/sample)');
ylabel('Phase ∠Y [rad]');
title('Phase Spectrum of y(n)');
saveas(fig2, fullfile(fig_dir, 'y_spectrum.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);
fprintf('Saved numeric data to: %s\n', data_dir);

function [fL, fR, BW, pkdB, fpk] = minus3dB_bw(f, mag, f0, search_width_bins)
% Finds -3 dB bandwidth around the local max nearest f0.
    if nargin < 4, search_width_bins = max(10, round(0.02*numel(f))); end
    magdB = 20*log10(mag + eps);
    [~, k0] = min(abs(f - f0));
    kL = max(1, k0 - search_width_bins);
    kR = min(numel(f), k0 + search_width_bins);
    [pkdB, kpk_rel] = max(magdB(kL:kR));
    kpk = kL + kpk_rel - 1;
    fpk = f(kpk);
    target = pkdB - 3;

    % left crossing
    k = kpk;
    while k > 1 && magdB(k) > target, k = k - 1; end
    if k == 1
        fL = f(1);
    else
        fL = interp1([magdB(k), magdB(k+1)], [f(k), f(k+1)], target, 'linear');
    end

    % right crossing
    k = kpk;
    while k < numel(f) && magdB(k) > target, k = k + 1; end
    if k == numel(f)
        fR = f(end);
    else
        fR = interp1([magdB(k-1), magdB(k)], [f(k-1), f(k)], target, 'linear');
    end

    BW = fR - fL;
end
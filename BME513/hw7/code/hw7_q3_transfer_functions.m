%% hw7_q3_transfer_functions.m
% TF estimates via cross-spectra

hw7_load_data;

x_list = {bp,     bp,      etco2,   etco2};
y_list = {fv,     hr,      fv,      hr};
labels = {'bp->fv','bp->hr','etco2->fv','etco2->hr'};

Lxy   = maxlag;
Nfft  = 2048;
faxis = (0:Nfft-1).' * (fs/Nfft);

Hhat_all = cell(1,4);

figure;
for k = 1:4
    x = x_list{k};
    y = y_list{k};
    
    % auto- and cross-correlations (biased)
    [rxx, ~] = xcorr(x, Lxy, 'biased');
    [ryx, ~] = xcorr(y, x, Lxy, 'biased');
    
    % spectra
    Sxx = fft(rxx, Nfft);
    Syx = fft(ryx, Nfft);
    
    % TF estimate
    Hhat = Syx ./ (Sxx + eps);
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
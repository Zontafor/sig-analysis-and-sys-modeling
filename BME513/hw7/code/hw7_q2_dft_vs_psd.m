%% hw7_q2_dft_vs_psd.m
% |DFT| vs PSD (via autocorrelation) for 4 datasets

hw7_load_data;

vars  = {bp, fv, hr, etco2};
names = {'bp', 'fv', 'hr', 'etco2'};

Lac   = maxlag;
Nfft  = 2048;
faxis = (0:Nfft-1).' * (fs/Nfft);

figure;
for k = 1:4
    x = vars{k};
    
    % raw DFT magnitude
    X    = fft(x, Nfft);
    magX = abs(X);
    
    % autocorrelation-based PSD (biased, max-lag Lac)
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
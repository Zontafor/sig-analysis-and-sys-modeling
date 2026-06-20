%% hw7_q4_irf_from_tf_and_ls.m
% IRFs via inverse FFT and LS FIR

hw7_load_data;
hw7_q3_transfer_functions;  % Hhat_all, x_list, y_list, labels

Lirf  = 40;                     % IRF length (samples)
Nfft  = 2048;
t_irf = (0:Lirf-1).' * Ts;

irf_ifft = cell(1,4);
irf_ls   = cell(1,4);

figure;
for k = 1:4
    x    = x_list{k};
    y    = y_list{k};
    Hhat = Hhat_all{k};
    
    % IFFT-based IRF (circular; take first Lirf taps)
    h_circ       = real(ifft(Hhat, Nfft));
    h_ifft       = h_circ(1:Lirf);
    irf_ifft{k}  = h_ifft;
    
    % LS FIR: y[n] ≈ sum_{m=0}^{Lirf-1} h[m] x[n-m]
    n0   = Lirf;
    nmax = numel(x)-1;
    M    = nmax - n0 + 1;
    Xreg = zeros(M, Lirf);
    yreg = zeros(M, 1);
    for i = 1:M
        n = n0 + i - 1;
        Xreg[i,:] = flipud(x(n-Lirf+1:n));
        yreg[i]   = y(n);
    end
    
    h_ls       = Xreg \ yreg;
    irf_ls{k}  = h_ls;
    
    subplot(2,2,k);
    plot(t_irf, h_ifft, 'LineWidth', 1.5); hold on;
    plot(t_irf, h_ls,   '--', 'LineWidth', 1.5);
    xlabel('time lag (s)');
    ylabel('h(t)');
    title(sprintf('IRF: %s', labels{k}));
    legend('IFFT(H)','LS-FIR','Location','best');
    grid on;
end
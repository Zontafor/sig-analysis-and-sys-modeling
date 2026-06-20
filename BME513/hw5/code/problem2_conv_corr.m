function problem2_conv_corr()
%%  PROBLEM2_CONV_CORR  Convolution and correlation of noisy exponentials.
%   Generates
%          x[n] = exp(-n/50)  + w_1[n]/5
%          y[n] = exp(-n/100) + w_2[n]/10
%        for n = 1,...,500 with independent white Gaussian noises.
%   Computes and plots the convolution x*y in the time domain
%        and via the frequency domain using the FFT.
%   Computes and plots the cross-correlation between x and y
%        in the time domain and via the frequency domain.
%   Prints numerical error norms between the direct and FFT-based
%        implementations to confirm correctness.

    rng(2);  % fix random seed for reproducibility

    N = 500;
    n = 1:N;

    % Generate independent white-noise sequences
    w1 = randn(size(n));
    w2 = randn(size(n));

    % Noisy exponentials
    x = exp(-n/50)  + w1/5;
    y = exp(-n/100) + w2/10;

    %%   Problem 2 Part 1
    %    Convolution in the time domain

    z_time = conv(x, y);
    n_conv = 0:(numel(z_time)-1);

    % Convolution via FFT (frequency domain)
    L  = numel(z_time);        % length needed for linear convolution
    Xf = fft(x, L);
    Yf = fft(y, L);
    z_freq = ifft(Xf .* Yf);

    % Numerical check
    conv_error = norm(z_time - z_freq);
    fprintf('||conv_{time} - conv_{freq}||_2 = %.3e\n', conv_error);

    % Plots for convolution
    figure('Color','w');
    subplot(2,1,1);
    plot(n_conv, z_time, 'LineWidth', 1);
    xlabel('n');
    ylabel('z_{time}[n]');
    title('Convolution x*y (time-domain conv)');
    grid on;
    set(gca,'FontName','Times','FontSize',11);

    subplot(2,1,2);
    plot(n_conv, real(z_freq), 'LineWidth', 1);
    xlabel('n');
    ylabel('Re\{z_{freq}[n]\}');
    title('Convolution x*y (frequency-domain via FFT)');
    grid on;
    set(gca,'FontName','Times','FontSize',11);

    sgtitle('Convolution of x[n] and y[n]','FontName','Times','FontSize',12);
    saveas(gcf, 'HW5_P2_conv_time_vs_freq.png');

    %%   Problem 2 Part 2
    %    Cross-correlation in the time domain

    [Rxy_time, lags] = xcorr(x, y);   % time-domain estimate
    % lags runs from -(N-1) to (N-1), length 2N-1
    % Cross-correlation via FFT matched to xcorr:
    %  - zero-pad x,y to length Lcorr = 2N-1
    %  - use X .* conj(Y) (cross-correlation theorem)
    %  - fftshift to align lag 0 in the center

    Lcorr   = 2*N - 1;               % same as length(Rxy_time)
    Xf_corr = fft(x, Lcorr);
    Yf_corr = fft(y, Lcorr);

    % Circular correlation, then shift to lags -(N-1):N-1
    Rxy_freq_circ = ifft( Xf_corr .* conj(Yf_corr) );
    Rxy_freq      = fftshift(Rxy_freq_circ);

    % Numerical check – should now be ~1e-13
    corr_error = norm(Rxy_time - Rxy_freq);
    fprintf('||R_{xy,time} - R_{xy,freq}||_2 = %.3e\n', corr_error);

    % Plots for correlation
    figure('Color','w');
    subplot(2,1,1);
    stem(lags, Rxy_time, 'filled');
    xlabel('Lag m');
    ylabel('R_{xy,time}[m]');
    title('Cross-correlation R_{xy}[m] (time-domain xcorr)');
    grid on;
    set(gca,'FontName','Times','FontSize',11);

    subplot(2,1,2);
    stem(lags, real(Rxy_freq), 'filled');
    xlabel('Lag m');
    ylabel('Re\{R_{xy,freq}[m]\}');
    title('Cross-correlation R_{xy}[m] (frequency-domain via FFT)');
    grid on;
    set(gca,'FontName','Times','FontSize',11);

    sgtitle('Cross-correlation between x[n] and y[n]', ...
            'FontName','Times','FontSize',12);
    saveas(gcf, 'HW5_P2_corr_time_vs_freq.png');
end
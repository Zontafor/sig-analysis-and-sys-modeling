function problem1_signal_analysis()
%%  Signal generation and spectral analysis on noisy sinusoids
%   Generates x[n] = cos(2*pi*n/250) + sin(2*pi*n/100 + pi/4) + c*w[n]
%        for n = 1,...,1000 with tunable noise level c.
%   Plots the time-domain signal.
%   Computes the N-point DFT and plots real, imaginary,
%        magnitude, and phase.
%   Computes and plots the biased autocorrelation estimate.
%   Sweeps different values of c to illustrate when periodicities
%        are visually discernible in time and frequency.

    rng(1);  % fix random seed for reproducibility

    % Basic parameters
    N  = 1000;
    n  = 1:N;
    fs = 1;          % normalized sampling frequency (arbitrary units)

    %%   Problem 1 Part 1

    c  = 1;
    w  = randn(size(n));  % white Gaussian noise, mean 0, variance 1
    x  = cos(2*pi*n/250) + sin(2*pi*n/100 + pi/4) + c*w;

    % Time-domain plot
    figure('Color','w');
    plot(n, x, 'LineWidth', 1);
    xlabel('n');
    ylabel('x[n]');
    title('Time-domain signal x[n] for c = 1');
    grid on;
    set(gca,'FontName','Times','FontSize',11);
    saveas(gcf, 'HW5_P1_time_series_c1.png');

    %%   Problem 1 Part 2
    %    DFT: real/imag, magnitude, and phase

    X = fft(x);
    k = 0:N-1;                 % DFT index
    f = k / N * fs;            % normalized frequency axis

    figure('Color','w');
    subplot(2,2,1);
    stem(f, real(X), 'filled'); grid on;
    xlabel('Normalized frequency');
    ylabel('Re\\{X[k]\\}');
    title('Real\\{X[k]\\}');
    set(gca,'FontName','Times','FontSize',10);

    subplot(2,2,2);
    stem(f, imag(X), 'filled'); grid on;
    xlabel('Normalized frequency');
    ylabel('Im\\{X[k]\\}');
    title('Imag\\{X[k]\\}');
    set(gca,'FontName','Times','FontSize',10);

    subplot(2,2,3);
    stem(f, abs(X), 'filled'); grid on;
    xlabel('Normalized frequency');
    ylabel('|X[k]|');
    title('Magnitude spectrum');
    set(gca,'FontName','Times','FontSize',10);

    subplot(2,2,4);
    stem(f, angle(X), 'filled'); grid on;
    xlabel('Normalized frequency');
    ylabel('\\angle X[k] (rad)');
    title('Phase spectrum');
    set(gca,'FontName','Times','FontSize',10);

    sgtitle('DFT of x[n] for c = 1','FontName','Times','FontSize',12);
    saveas(gcf, 'HW5_P1_DFT_components_c1.png');

    %%   Problem 1 Part 3
    %    Autocorrelation estimate

    [rxx, lags] = xcorr(x, 'biased');

    figure('Color','w');
    stem(lags, rxx, 'filled');
    xlabel('Lag m');
    ylabel('\\hat{R}_{xx}[m]');
    title('Biased autocorrelation estimate of x[n] (c = 1)');
    grid on;
    set(gca,'FontName','Times','FontSize',11);
    saveas(gcf, 'HW5_P1_autocorr_c1.png');

    %%   Problem 1 Part 3
    %    Sweep different values of c to illustrate visibility of periodicities in time and frequency domains

    c_vals = [0.25, 0.5, 1, 2, 4];
    numC   = numel(c_vals);

    % Time-domain sweep
    figure('Color','w');
    for idx = 1:numC
        c = c_vals(idx);
        w = randn(size(n));
        x = cos(2*pi*n/250) + sin(2*pi*n/100 + pi/4) + c*w;
        subplot(numC,1,idx);
        plot(n, x, 'LineWidth', 1);
        ylabel(sprintf('c = %.2f', c));
        grid on;
        set(gca,'FontName','Times','FontSize',9);
        if idx == 1
            title('Time-domain x[n] for different noise levels c');
        end
        if idx == numC
            xlabel('n');
        end
    end
    saveas(gcf, 'HW5_P1_time_sweep_c_vals.png');

    % Frequency-domain sweep (magnitude only, low-frequency zoom)
    figure('Color','w');
    for idx = 1:numC
        c = c_vals(idx);
        w = randn(size(n));
        x = cos(2*pi*n/250) + sin(2*pi*n/100 + pi/4) + c*w;
        X = fft(x);
        magX = abs(X);

        subplot(numC,1,idx);
        stem(f, magX, 'filled');
        xlim([0 0.2]); % zoom near DC to highlight sinusoidal lines
        grid on;
        ylabel(sprintf('c = %.2f', c));
        set(gca,'FontName','Times','FontSize',9);
        if idx == 1
            title('Magnitude |X[k]| for different noise levels c');
        end
        if idx == numC
            xlabel('Normalized frequency');
        end
    end
    saveas(gcf, 'HW5_P1_mag_sweep_c_vals.png');
end
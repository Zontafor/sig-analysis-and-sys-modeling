%% hw7_q1_hist_scatter.m
% histograms with Gaussian fits + bp vs fv scatter

hw7_load_data;

vars  = {bp, fv, hr, etco2};
names = {'bp', 'fv', 'hr', 'etco2'};
nbins = 20;

figure;
for k = 1:4
    x = vars{k};
    subplot(2,2,k);
    
    % histogram (PDF normalized)
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

% bp vs fv scatter
figure;
plot(bp, fv, '.', 'MarkerSize', 6);
xlabel('bp (mmHg, demeaned)');
ylabel('fv (cm/s, demeaned)');
title('Contemporaneous bp vs fv');
grid on; axis tight;
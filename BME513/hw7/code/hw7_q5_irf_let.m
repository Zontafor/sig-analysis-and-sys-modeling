%% hw7_q5_irf_let.m
% IRFs via Laguerre Expansion Technique (LET) for 4 pairs

hw7_load_data;

x_list = {bp,     bp,      etco2,   etco2};
y_list = {fv,     hr,      fv,      hr};
labels = {'bp->fv','bp->hr','etco2->fv','etco2->hr'};

Llag   = 4;       % number of Laguerre functions (per original code)
Mlag   = 50;      % memory length for Laguerre expansion
alpha  = 0.5;     % alpha parameter (a = sqrt(alpha))
Lirf   = Mlag;    % IRF length in samples (same as M)
t_irf  = (0:Lirf-1).' * Ts;

h_let  = cell(1,4);
Q_let  = zeros(4,1);

figure;
for k = 1:4
    x = x_list{k};
    y = y_list{k};
    
    [hest, yest, Gest, Q, LC] = LET1_single_pair(x, y, Llag, Mlag, alpha);
    h_let{k} = hest;
    Q_let(k) = Q;
    
    subplot(2,2,k);
    plot(t_irf, hest, 'LineWidth', 1.5);
    xlabel('time lag (s)');
    ylabel('h_{LET}(t)');
    title(sprintf('LET IRF: %s (Q=%.3f)', labels{k}, Q));
    grid on;
end
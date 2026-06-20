%% hw7_q4_diagnostics.m

pairs = labels;  % {'bp->fv','bp->hr','etco2->fv','etco2->hr'}

fprintf('\nQ4 Diagnostics\n');

for k = 1:4
    x = x_list{k};
    y = y_list{k};

    % rebuild LS regression matrix for condition number
    n0   = Lirf;
    nmax = numel(x)-1;
    Mls  = nmax - n0 + 1;
    Xreg = zeros(Mls, Lirf);
    for i = 1:Mls
        n = n0 + i - 1;
        Xreg(i,:) = flipud(x(n-Lirf+1:n));
    end

    condX = cond(Xreg' * Xreg);

    % energy and roughness of IRFs
    h_ls   = irf_ls{k};
    h_ifft = irf_ifft{k};
    h_letk = h_let{k}(1:Lirf);          % match length

    E_ls   = sum(h_ls.^2);
    E_ifft = sum(h_ifft.^2);
    E_let  = sum(h_letk.^2);

    % second-difference roughness
    R_ls   = sum(diff(h_ls,2).^2);
    R_ifft = sum(diff(h_ifft,2).^2);
    R_let  = sum(diff(h_letk,2).^2);

    fprintf('\nPair %s:\n', pairs{k});
    fprintf('  cond(X''X)     = %.2e\n', condX);
    fprintf('  Energy ratio   E_LS / E_LET   = %.2f\n', E_ls / max(E_let,eps));
    fprintf('  Roughness R_LS / R_LET        = %.2f\n', R_ls / max(R_let,eps));

    if condX > 1e6 || E_ls/E_let > 10 || R_ls/R_let > 10
        fprintf('  -> WARNING: LS FIR IRF likely ill-conditioned / overfitted.\n');
    else
        fprintf('  -> LS FIR IRF appears numerically well-conditioned.\n');
    end
end
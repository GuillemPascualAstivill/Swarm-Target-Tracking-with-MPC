function exact_W2 = compute_continuous_W2(agents, grid_res)
    % Parameters for "Exact" Evaluation
    u = 1e-5;         % Very low regularization to mimic unregularized OT
    realmin = 1e-200;
    N = 1000;
    c= 0.05;
    h = c * N^(-1/(2+2));
    ST_eval = 2000;
    [X, Y] = meshgrid(linspace(0, 1, grid_res), linspace(0, 1, grid_res));
    dx = 1 / (grid_res - 1);
    sigma_kde_pixels = h / dx;
    sigma_sink_pixels = sqrt(u/ 2) / dx;
    edges = linspace(0, 1, grid_res + 1);
    %Target Definition (change for target being used)
    theta = atan2(Y - 0.5, X - 0.5);
    r_grid = sqrt((X - 0.5).^2 + (Y - 0.5).^2);
    r_star = 0.25 + 0.1 * cos(5 * theta);
    target_thickness = 0.04;
    target_distr = exp(-(r_grid - r_star).^2 / (2 * target_thickness^2)) + realmin;
    target_distr = target_distr / sum(target_distr(:));
    final_counts = histcounts2(agents(:,2), agents(:,1), edges, edges);
    KDE = imgaussfilt(final_counts, sigma_kde_pixels, 'Padding', 'symmetric');
    KDE = KDE / sum(KDE(:));


    % --- Sinkhorn Setup ---
    % Kernel bandwidth for the evaluation
    sigma_sink_pixels = sqrt(u / 2) / dx;
    K_blur = @(M) imgaussfilt(M, sigma_sink_pixels, 'Padding', 'symmetric');

    % Initialize dual scaling factors
    b_sink = ones(grid_res, grid_res);
    a_sink = ones(grid_res, grid_res);

    % --- Alternating Projections ---
    for i = 1:ST_eval
        a_sink = KDE ./ max(K_blur(b_sink), realmin);
        b_sink = target_distr ./ max(K_blur(a_sink), realmin);
    end

    % --- Compute Wasserstein-2 Distance ---
    phi = u * log(max(a_sink, realmin));
    psi = u * log(max(b_sink, realmin));

    % Summing over the grid (Approximating the integral \int phi d_rho)
    W2_squared = sum(phi(:) .* KDE(:)) + sum(psi(:) .* target_distr(:));

    % Correct for potential numerical offsets and return W2
    exact_W2 = sqrt(max(W2_squared, 0));
end

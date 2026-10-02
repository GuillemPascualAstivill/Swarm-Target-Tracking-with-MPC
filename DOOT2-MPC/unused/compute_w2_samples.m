function W2 = compute_w2_samples(final_samples, target_samples)
% Computes empirical Wasserstein-2 distance between two sets of 2D samples
%
% Inputs:
%   final_samples  : N x 2 matrix
%   target_samples : N x 2 matrix
%
% Output:
%   W2 : Wasserstein-2 distance

    % Dimension checks
    if size(final_samples,2) ~= 2 || size(target_samples,2) ~= 2
        error('Samples must be N x 2 matrices');
    end

    if size(final_samples,1) ~= size(target_samples,1)
        error('Both sets must have the same number of samples');
    end

    N = size(final_samples,1);

    % Squared Euclidean cost matrix
    C = pdist2(final_samples, target_samples).^2;

    % Hungarian assignment
    pairs = matchpairs(C, 1e10);

    % Compute transport cost manually
    total_cost = 0;
    for k = 1:size(pairs,1)
        i = pairs(k,1);
        j = pairs(k,2);
        total_cost = total_cost + C(i,j);
    end

    % Wasserstein distance
    W2 = sqrt(total_cost / N);

end

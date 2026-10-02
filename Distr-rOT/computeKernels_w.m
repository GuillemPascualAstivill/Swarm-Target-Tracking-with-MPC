function kernels = computeKernels_w(pos, params)
%computeKernels_w  Agent-agent Wendland kernels.
%   kernels.K_agent   density kernel (h_agent):  rho_N(x_i) = sum_j K_agent(i,j) / N
%   kernels.K_field   potential basis (h_field): phi(x_i)   = K_field(i,:) * alpha
%   kernels.D_x, D_y  basis gradients:           grad phi(x_i) = [D_x(i,:), D_y(i,:)] * alpha

diffs = permute(pos, [1 3 2]) - permute(pos, [3 1 2]);   % diffs(i,j,:) = x_i - x_j
sq_dists = sum(diffs.^2, 3);

kernels.K_agent = wendland_kernel(sq_dists, params.h_agent);
[kernels.K_field, Dfac] = wendland_kernel(sq_dists, params.h_field);
kernels.D_x = Dfac .* diffs(:, :, 1);
kernels.D_y = Dfac .* diffs(:, :, 2);

end

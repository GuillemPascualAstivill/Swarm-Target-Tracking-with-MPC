%% check_ridge.m -- does the planner converge? One planning step (initial
% swarm, target at phase 0), with and without the ridge, at 3000 and 30000
% primal-dual iterations. Run parameters_w.m first. Takes about a minute.

pos = params.pos_cont0;
P = params;
P.rho_target = star_target_w(params, 0);
kernels  = computeKernels_w(pos, P.ant_cords, P);
integral = computeTargetint_w(pos, P.ant_cords, P);
N = P.N;
w = P.w_gain * (sum(integral.agents_ag(1:N, 1:N), 2) / N - integral.target_ag(1:N));
c_sq = conformal(pos, P.ant_cords, P);
velocity = @(a) -[kernels.D_x(1:N, :) * a, kernels.D_y(1:N, :) * a];

for d = [0, 4e-3]
    P.delta_ridge = d;
    P.pd_iters = 3000;   v1 = velocity(primal_dual(kernels, integral, w, c_sq, P));
    P.pd_iters = 30000;  v2 = velocity(primal_dual(kernels, integral, w, c_sq, P));
    fprintf('delta_ridge = %g: 3000 vs 30000 iterations differ by %.1f%%, max speed %.3f\n', ...
        d, 100 * norm(v1 - v2, 'fro') / norm(v2, 'fro'), max(vecnorm(v1, 2, 2)));
end

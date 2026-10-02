function integral = computeTargetint_w(pos, rho_target, params)
%computeTargetint_w  Grid integrals coupling the agents to the target.
%   integral.agents_ag(j,i)          = int K_field(x - x_j) K_agent(x - x_i) dx
%   integral.target_ag(j)            = int K_field(x - x_j) rho_target(x) dx
%   integral.rho_target_at_agents(i) = rho_target(x_i)                          ('kl')
%   integral.K_reg(i,j)              = int K_field(x - x_i) K_field(x - x_j) dx  ('l2')

dx = params.dx;
dy = params.dy;

diffs = permute(params.grid_pts, [1 3 2]) - permute(pos, [3 1 2]);
sq_dists = sum(diffs.^2, 3);                             % grid points x agents
K_grid_agent = wendland_kernel(sq_dists, params.h_agent);
K_grid_field = wendland_kernel(sq_dists, params.h_field);

integral.target_ag = K_grid_field' * (rho_target * dx * dy);
integral.agents_ag = (K_grid_field' * K_grid_agent) * dx * dy;

% Terms used only by the chosen terminal cost
switch params.terminal
    case 'kl'
        integral.rho_target_at_agents = interp2(params.X, params.Y, ...
            reshape(rho_target, size(params.X)), pos(:,1), pos(:,2), 'linear', 0);
    case 'l2'
        integral.K_reg = (K_grid_field' * K_grid_field) * dx * dy;
    otherwise
        error('params.terminal must be ''kl'' or ''l2''.');
end

end

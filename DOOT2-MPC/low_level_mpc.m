function [u, n_active] = low_level_mpc(pos, vel, p_ref, v_tgt, d, mpc)
%LOW_LEVEL_MPC  One receding-horizon step of the low-level MPC for all agents.
%   Rows are agents and columns are axes; every agent-axis is an independent
%   problem for the plant xdot = v + d, vdot = u (see LOW_LEVEL_MPC_SETUP).
%   Each agent tracks the point p_ref moving at v_tgt, with error state
%       zeta = [pos - p_ref; vel - v_ref],   v_ref = v_tgt - d
%   so zeta = 0 is an equilibrium of the plant in the wind d. Here d is the
%   wind the agent believes in (measured, estimated or zero); if the true
%   wind differs by a constant, the gain leaves a standing offset
%   (k_v/k_p) (d_true - d).
%
%   u (N x 2) is the applied acceleration: u = -K zeta if mpc.u_max is
%   empty, otherwise the first move of the box-constrained QP over the
%   horizon. n_active counts the agent-axes whose QP has an active bound.

    e_p = pos - p_ref;
    e_v = vel - (v_tgt - d);

    if isempty(mpc.u_max)
        u = -mpc.K(1) * e_p - mpc.K(2) * e_v;
        n_active = 0;
    else
        % One QP per agent-axis (column). Where the unconstrained optimum
        % -H\f satisfies the bounds over the whole horizon it is also the QP
        % solution, so quadprog only runs where a bound is active.
        f = mpc.F * [e_p(:)'; e_v(:)'];
        U = -(mpc.H \ f);
        active = find(any(abs(U) > mpc.u_max, 1));
        for j = active
            [U(:, j), ~, flag] = quadprog(mpc.H, f(:, j), [], [], [], [], mpc.lb, mpc.ub, ...
                min(max(U(:, j), mpc.lb), mpc.ub), mpc.qp_opts);
            if flag <= 0
                error('low_level_mpc:qp', 'quadprog failed with exitflag %d.', flag);
            end
        end
        u = reshape(U(1, :), size(pos));
        n_active = numel(active);
    end
end

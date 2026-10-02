function [u_opt, v_next, x_next, info] = low_level_mpc(x_m, v_m, x_nom, v_tgt, t, wind_fn, mpc, U_ws, d_hat)
% LOW_LEVEL_MPC  One receding-horizon step of the per-agent low-level MPC.
%
%   Drives the agent toward the destination x_nom under second-order
%   dynamics with drift (xdot = v + d, vdot = u) and returns the
%   dynamically consistent triple (x_next, v_next, u_opt). One step of size
%   mpc.dt. Called per agent, per inner step.
%
%   Error state -- both rows in the GROUND frame:
%
%     zeta = [ p - x_nom ;  v - v_ref ],    v_ref = v_tgt - d
%
%   so the second row is (v + d) - v_tgt = agent ground velocity minus the
%   velocity of the target point. Equivalently zeta = [e; edot] with
%   e = p - x_nom(t). The -d is not a heuristic feedforward: it makes
%   zeta = 0 a genuine equilibrium of the plant. Without it, a persistent
%   drift leaves a standing offset e_inf = (k_v/k_p)*d = sqrt(2)*d/omega,
%   which grows as the loop is tuned smoother.
%
%   v_tgt is the velocity OF THE TARGET POINT, and the caller owns it:
%     v_tgt = [0 0]    x_nom frozen over the inner loop (setpoint regulator)
%     v_tgt = v_nom    x_nom advancing at the transport velocity (ramp)
%   Note v_nom is DOOT's instantaneous velocity command, NOT a required
%   arrival velocity, which is why it appears here only through v_tgt and
%   carries no weight in the cost.
%
%   Three paths, selected by mpc.constrained and mpc.mode:
%     unconstrained    : u = -K*zeta, the constant gain from
%                        LOW_LEVEL_MPC_SETUP. This is LQR.
%     'clamp'          : the same gain, saturated after the fact. Cheap,
%                        but saturated LQR is not the constrained optimum.
%     'qp'             : solve the condensed box-constrained QP over the
%                        horizon and apply its first move. The solution is
%                        piecewise affine rather than linear, which is what
%                        makes this an MPC.
%   All three coincide while the box is slack, so a null result in the
%   clamp-vs-QP comparison may mean the constraint never activated.
%
%   Inputs
%     x_m      : current position (1x2)
%     v_m      : current velocity (1x2, persistent state)
%     x_nom    : destination position (1x2); held fixed across the inner loop
%     v_tgt    : velocity of the target point (1x2); [0 0] for a frozen x_nom
%     t        : current time (passed to the drift query)
%     wind_fn  : drift handle, d = wind_fn(t, x) (world frame; planar
%                components used)
%     mpc      : struct from LOW_LEVEL_MPC_SETUP
%     U_ws     : warm start for the QP path, Nh x 2 (optional). Pass back
%                info.U_next from the previous inner step; in a receding
%                horizon the solution barely moves, so this cuts iterations
%                sharply. Ignored unless mode is 'qp'.
%
%   Outputs
%     u_opt  : first applied acceleration (1x2)
%     v_next : next velocity  = v_m + dt*u_opt        (1x2)
%     x_next : next position  = x_m + dt*(v_next + d) (1x2, true drift)
%     info   : (optional) x_nom, v_ref, zeta0, drift d, u_opt, and for the
%              QP path U (the solved sequence), U_next (shifted, to feed
%              back as U_ws) and n_active (how many moves sit on the box)

    dt = mpc.dt;
    kp = mpc.K(1);
    kv = mpc.K(2);

    % Row vectors so the per-axis operations broadcast cleanly
    x_m   = x_m(:).';
    v_m   = v_m(:).';
    x_nom = x_nom(:).';

    if nargin < 4 || isempty(v_tgt)
        v_tgt = [0, 0];
    else
        v_tgt = v_tgt(:).';
    end

    % Drift at the current position (re-queried each step); use planar part
    d = wind_fn(t, x_m);
    d = d(:).';
    d = d(1:2);

    % What the controller believes the drift to be. The true d is always
    % returned in info for the plant integration; only the control law is
    % affected here.
    %   'measured'  - the agent reads the true local wind (assumes an
    %                 anemometer); zeta = 0 is exactly an equilibrium.
    %   'none'      - the agent knows nothing of d and sees it only through
    %                 the state it measures next step.
    %   'estimated' - the agent uses d_hat, inferred by the caller from the
    %                 gap between predicted and measured motion. Needs only
    %                 self-localization, not a wind sensor.
    if isfield(mpc, 'drift_mode')
        drift_mode = mpc.drift_mode;
    elseif isfield(mpc, 'drift_feedforward') && ~mpc.drift_feedforward
        drift_mode = 'none';
    else
        drift_mode = 'measured';
    end

    switch drift_mode
        case 'measured'
            d_ctrl = d;
        case 'none'
            d_ctrl = [0, 0];
        case 'estimated'
            if nargin < 9 || isempty(d_hat)
                d_ctrl = [0, 0];
            else
                d_hat  = d_hat(:).';
                d_ctrl = d_hat(1:2);
            end
        otherwise
            error('low_level_mpc:driftMode', ...
                  'drift_mode must be measured, none or estimated.');
    end

    v_ref = v_tgt - d_ctrl;
    e_p   = x_m - x_nom;
    e_v   = v_m - v_ref;

    constrained = isfield(mpc, 'constrained') && mpc.constrained;
    use_qp_path = constrained && isfield(mpc, 'mode') && strcmp(mpc.mode, 'qp');

    U_opt    = [];
    U_next   = [];
    n_active = 0;

    if use_qp_path
        Nh = mpc.Nh;
        if nargin < 8 || isempty(U_ws)
            U_ws = zeros(Nh, 2);
        end
        U_opt = zeros(Nh, 2);
        u_opt = zeros(1, 2);
        for a = 1:2
            % Only the linear term depends on the state; H is precomputed.
            f = mpc.qp.F * [e_p(a); e_v(a)];
            U = quadprog(mpc.qp.H, f, [], [], [], [], ...
                         mpc.qp.lb, mpc.qp.ub, U_ws(:, a), mpc.qp.opts);
            if isempty(U) || any(~isfinite(U))
                % Solver failed. Fall back to the clamped gain rather than
                % propagating a NaN through the swarm, and say so.
                warning('low_level_mpc:qpFallback', ...
                        'quadprog failed on axis %d; using clamped gain.', a);
                U = max(min(-kp*e_p(a) - kv*e_v(a), mpc.u_max), -mpc.u_max) ...
                    * ones(Nh, 1);
            end
            U_opt(:, a) = U;
            u_opt(a)    = U(1);
        end
        % Shift for the next step's warm start, repeating the last move.
        U_next   = [U_opt(2:end, :); U_opt(end, :)];
        n_active = sum(abs(U_opt(:)) >= mpc.u_max - 1e-12);
    else
        u_opt = -kp * e_p - kv * e_v;
        if constrained
            % Saturated LQR, not the constrained optimum. See the header.
            u_opt = max(min(u_opt, mpc.u_max), -mpc.u_max);
            n_active = sum(abs(u_opt) >= mpc.u_max - 1e-12);
        end
    end

    % Integrate: velocity, then ground position including the true drift
    v_next = v_m + dt * u_opt;
    x_next = x_m + dt * (v_next + d);

    if nargout > 3
        info = struct('x_nom', x_nom, ...
                      'v_ref', v_ref, ...
                      'zeta0', [e_p; e_v], ...
                      'd',     d, ...
                      'u_opt', u_opt, ...
                      'U',        U_opt, ...
                      'U_next',   U_next, ...
                      'n_active', n_active);
    end
end

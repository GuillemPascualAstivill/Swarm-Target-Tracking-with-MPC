function mpc = low_level_mpc_setup(dt, Q_x, Q_u, Nh, gamma, u_max, Q_v, mode)
% LOW_LEVEL_MPC_SETUP  Precompute the low-level MPC gain and terminal cost.
%
%   Finite-horizon LQR on the per-axis 2-state double integrator (state
%   [p; v], control u = acceleration). The stage cost tracks position, and
%   optionally the velocity error too (Q_v). The terminal cost is
%   M = gamma * P*, where P* is the infinite-horizon cost-to-go from
%   idare. Called once, before the outer DOOT loop.
%
%   Cost (per axis, zeta = [p - x_nom; v - v_ref]):
%
%     J = sum_{k=0}^{Nh-1} [ Q_x*(p_k - x_nom)^2 + Q_v*(v_k - v_ref)^2
%                            + Q_u*u_k^2 ]
%         + zeta_Nh' * M * zeta_Nh,      M = gamma * P*
%
%   Why Q_v: with Q_v = 0 the closed loop has damping ratio 1/sqrt(2)
%   whatever Q_x/Q_u is, so the only knob is the bandwidth
%   omega = (Q_x/Q_u)^(1/4). A nonzero Q_v decouples damping from
%   bandwidth, and prices the velocity error directly -- which is the
%   quantity the macroscopic bound charges for. It also breaks the closed
%   form for P and K that the paper's appendix quotes. Default 0, which
%   reproduces the position-only behavior exactly.
%
%   Why gamma: M = P* (gamma = 1) is the exact cost-to-go, so the backward
%   recursion sits at its fixed point, the gain is the constant idare gain
%   and Nh is inert. gamma > 1 overprices the terminal state, which raises
%   the gain, speeds arrival and shrinks the pace deficit, at the cost of
%   smoothness and some overshoot. gamma >= 1 keeps the CLF condition
%   R(M) <= M, hence the standard MPC stability argument (Limon et al.,
%   Automatica 2006: a scaled terminal cost enlarges the domain of
%   attraction and can replace the terminal constraint set).
%
%   Inputs
%     dt    : inner discretization step            (sets A, B)
%     Q_x   : position-tracking weight             (Q(1,1))
%     Q_u   : control-effort weight                (R)
%     Nh    : prediction horizon (optional, default 10). Load-bearing when
%             gamma ~= 1; inert when gamma = 1.
%     gamma : terminal-cost scale, M = gamma*P* (optional, default 1)
%     u_max : actuator box |u| <= u_max (optional). [] = unconstrained.
%     Q_v   : velocity-tracking weight, Q(2,2) (optional, default 0)
%     mode  : how u_max is enforced (optional, default 'clamp'):
%               'clamp' - saturate the unconstrained gain after the fact.
%                         Cheap, but it is saturated LQR, NOT the solution
%                         of the constrained problem, and the two differ
%                         precisely when the constraint is active.
%               'qp'    - solve the condensed box-constrained QP over the
%                         horizon and apply its first move. This is what
%                         makes the low level an MPC rather than an LQR.
%             Ignored when u_max is empty.
%
%   Output (struct mpc)
%     mpc.K      : feedback gain [k_p, k_v] (1x2) at k = 0, law u = -K*zeta
%     mpc.M      : terminal-cost matrix (2x2) = gamma * P*
%     mpc.P_inf  : idare solution P* (2x2), the gamma = 1 reference
%     mpc.P0     : P_0 from the backward recursion (optimal cost-to-go at k=0)
%     mpc.K_inf  : constant idare gain, for comparison against mpc.K
%     mpc.A,B    : discrete per-axis model
%     mpc.Q,R    : weights, Q = diag(Q_x, Q_v), R = Q_u
%     mpc.Nh,dt  : horizon and step (passed through)
%     mpc.gamma  : terminal-cost scale
%     mpc.u_max  : actuator limit
%     mpc.use_qp : true if u_max set (phase-2 QP path)

    if nargin < 4 || isempty(Nh),    Nh    = 10;  end
    if nargin < 5 || isempty(gamma), gamma = 1;   end
    if nargin < 6,                   u_max = []; end
    if nargin < 7 || isempty(Q_v),   Q_v   = 0;   end
    if nargin < 8 || isempty(mode),  mode  = 'clamp'; end

    mode = lower(mode);
    if ~ismember(mode, {'clamp', 'qp'})
        error('low_level_mpc_setup:mode', ...
              'mode must be ''clamp'' or ''qp'', got ''%s''.', mode);
    end

    if Q_v < 0
        error('low_level_mpc_setup:Qv', 'Q_v must be nonnegative.');
    end

    if gamma < 1
        warning('low_level_mpc_setup:gamma', ...
                ['gamma = %.3g < 1: M is cheaper than the true tail cost, ' ...
                 'so the CLF condition R(M) <= M is not guaranteed.'], gamma);
    end

    % Per-axis double integrator, ZOH discretized at dt
    A = [1, dt; 0, 1];
    B = [dt^2/2; dt];
    Q = diag([Q_x, Q_v]);       % Q_v = 0 recovers the position-only cost
    R = Q_u;

    % Infinite-horizon cost-to-go P* and its constant gain (the gamma = 1
    % reference). (A, Q^(1/2)) stays observable through the position channel
    % even with zero velocity weight, so P* exists and is positive definite.
    [P_inf, K_inf] = idare(A, B, Q, R);

    % Terminal cost, then the finite-horizon backward Riccati recursion:
    %   P_Nh = M
    %   K_k  = (R + B'P_{k+1}B) \ (B'P_{k+1}A)
    %   P_k  = Q + A'P_{k+1}A - A'P_{k+1}B*K_k
    % Only K_0 is applied (receding horizon), so the applied gain is still a
    % CONSTANT -- gamma and Nh change which constant, not its time-variance.
    M = gamma * P_inf;
    P = M;
    Kk = K_inf;                 % defensive default if Nh = 0
    for k = Nh-1:-1:0
        Kk = (R + B'*P*B) \ (B'*P*A);
        P  = Q + A'*P*A - A'*P*B*Kk;
    end

    % Phase-1 closed loop is LTI, so stability is checkable directly. (In
    % phase 2 the QP makes the law piecewise-affine and this check no longer
    % applies -- that is where the CLF/terminal-cost argument earns its keep.)
    rho = max(abs(eig(A - B*Kk)));
    if rho >= 1
        warning('low_level_mpc_setup:unstable', ...
                'Closed loop not Schur stable: max|eig(A-BK)| = %.4f.', rho);
    end

    mpc = struct();
    mpc.K      = Kk;         % 1x2, K_0
    mpc.M      = M;          % 2x2 terminal cost
    mpc.P_inf  = P_inf;      % 2x2
    mpc.P0     = P;          % 2x2
    mpc.K_inf  = K_inf;      % 1x2
    mpc.A      = A;
    mpc.B      = B;
    mpc.Q      = Q;
    mpc.R      = R;
    mpc.Q_x    = Q_x;
    mpc.Q_v    = Q_v;
    mpc.Nh     = Nh;
    mpc.dt     = dt;
    mpc.gamma  = gamma;
    mpc.rho_cl = rho;
    mpc.u_max  = u_max;
    % What the control law believes the drift to be: 'measured' (reads the
    % true local wind, the original behaviour), 'none' (knows nothing of it)
    % or 'estimated' (uses a d_hat the caller infers from the gap between
    % predicted and measured motion). See LOW_LEVEL_MPC.
    mpc.drift_mode = 'measured';
    mpc.drift_feedforward = true;   % deprecated; kept for older callers
    mpc.constrained = ~isempty(u_max);
    mpc.mode   = mode;
    mpc.use_qp = mpc.constrained;   % deprecated alias; use mpc.constrained

    % Condense the horizon problem once. H does not depend on the state, so
    % every agent, axis and inner step reuses it; only the linear term moves.
    if mpc.constrained && strcmp(mode, 'qp')
        mpc.qp = low_level_mpc_condense(A, B, Q, R, M, Nh, u_max);
    else
        mpc.qp = [];
    end

    % With gamma = 1 and no constraint the QP's first move IS -K*zeta, so
    % the two modes coincide until the box binds. Worth knowing when the
    % clamp-vs-QP comparison shows no difference: it may mean the
    % constraint never activated, not that the distinction does not matter.
    if mpc.constrained
        u_lqr_typ = norm(Kk);       % |u| for a unit-norm error state
        if u_max > 10 * u_lqr_typ
            warning('low_level_mpc_setup:slack', ...
                    ['u_max = %.3g is large next to |K| = %.3g; the box ' ...
                     'may never activate.'], u_max, u_lqr_typ);
        end
    end
end

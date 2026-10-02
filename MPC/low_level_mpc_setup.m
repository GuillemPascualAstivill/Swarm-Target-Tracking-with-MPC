function mpc = low_level_mpc_setup(dt, Q_x, Q_u, Nh, gamma, u_max)
%LOW_LEVEL_MPC_SETUP  Precompute the low-level MPC solved by LOW_LEVEL_MPC.
%   Per agent and axis: double integrator with step dt, error state
%   zeta = [p - p_ref; v - v_ref], control u = acceleration, and cost
%       J = sum_{k=0}^{Nh-1} ( Q_x (p_k - p_ref)^2 + Q_u u_k^2 ) + zeta_Nh' M zeta_Nh
%   with terminal cost M = gamma * P_inf, P_inf the infinite-horizon
%   cost-to-go from idare. gamma = 1 makes M the exact tail cost, so the gain
%   is the idare gain and Nh is inert; gamma > 1 raises the gain (faster
%   tracking, less smooth).
%
%   u_max = []   unconstrained (default). The optimum is the constant gain
%                u = -K zeta from the backward Riccati recursion.
%   u_max > 0    |u_k| <= u_max on each axis. Eliminating the states,
%                [zeta_1; ...; zeta_Nh] = Sx zeta_0 + Su U, leaves the QP in
%                U = [u_0; ...; u_{Nh-1}]
%                    min 0.5 U' H U + (F zeta_0)' U   s.t.  |U| <= u_max
%                which LOW_LEVEL_MPC solves at every step, applying u_0.

    if nargin < 4 || isempty(Nh),    Nh = 10;    end
    if nargin < 5 || isempty(gamma), gamma = 1;  end
    if nargin < 6,                   u_max = []; end
    if ~isempty(u_max) && ~(isscalar(u_max) && u_max > 0)
        error('u_max must be [] or a positive scalar.');
    end
    if gamma < 1
        warning('low_level_mpc_setup:gamma', ...
            'gamma = %.3g < 1: M underprices the tail cost, stability is not guaranteed.', gamma);
    end

    % Per-axis double integrator, ZOH discretized at dt
    A = [1, dt; 0, 1];
    B = [dt^2/2; dt];
    Q = diag([Q_x, 0]);
    R = Q_u;

    % Terminal cost, then the backward Riccati recursion from P_Nh = M. Only
    % K_0 is applied (receding horizon), so the unconstrained law is constant.
    [P_inf, K] = idare(A, B, Q, R);
    M = gamma * P_inf;
    P = M;
    for k = Nh-1:-1:0
        K = (R + B'*P*B) \ (B'*P*A);
        P = Q + A'*P*A - A'*P*B*K;
    end
    if max(abs(eig(A - B*K))) >= 1
        warning('low_level_mpc_setup:unstable', 'The unconstrained closed loop is not stable.');
    end

    mpc = struct('dt', dt, 'Nh', Nh, 'gamma', gamma, 'u_max', u_max, ...
        'A', A, 'B', B, 'Q', Q, 'R', R, 'M', M, 'K', K);

    if ~isempty(u_max)
        % Condense: zeta_k = A^k zeta_0 + sum_{j<k} A^(k-1-j) B u_j, k = 1..Nh
        Sx = zeros(2*Nh, 2);
        Su = zeros(2*Nh, Nh);
        for k = 1:Nh
            Sx(2*k-1:2*k, :) = A^k;
            for j = 1:k
                Su(2*k-1:2*k, j) = A^(k-j) * B;
            end
        end
        Qbar = blkdiag(kron(eye(Nh-1), Q), M);   % stage cost on zeta_1..zeta_{Nh-1}, M on zeta_Nh
        H = 2 * (Su' * Qbar * Su + R * eye(Nh));
        mpc.H = (H + H') / 2;
        mpc.F = 2 * Su' * Qbar * Sx;
        mpc.lb = -u_max * ones(Nh, 1);
        mpc.ub =  u_max * ones(Nh, 1);
        mpc.qp_opts = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
    end
end

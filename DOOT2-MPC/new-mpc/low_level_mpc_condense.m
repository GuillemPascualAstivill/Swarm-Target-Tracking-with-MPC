function qp = low_level_mpc_condense(A, B, Q, R, M, Nh, u_max)
% LOW_LEVEL_MPC_CONDENSE  Condensed box-constrained QP for the low-level MPC.
%
%   Eliminates the states from the finite-horizon problem so the only
%   decision variable is the input sequence U = [u_0; ...; u_{Nh-1}], one
%   axis at a time. Called once from LOW_LEVEL_MPC_SETUP; the result is
%   reused by every agent, every axis and every inner step, because the
%   Hessian does not depend on the state.
%
%   Problem (per axis, zeta = [p - x_nom; v - v_ref]):
%
%     min_U  sum_{k=0}^{Nh-1} zeta_k' Q zeta_k + R u_k^2 + zeta_Nh' M zeta_Nh
%     s.t.   zeta_{k+1} = A zeta_k + B u_k,   |u_k| <= u_max
%
%   The zeta_0 term is constant in U and is dropped. Substituting
%
%     zeta_k = A^k zeta_0 + sum_{j<k} A^{k-1-j} B u_j
%
%   for k = 1..Nh gives zeta = Sx*zeta_0 + Su*U and hence
%
%     J(U) = 0.5 U' H U + (F zeta_0)' U + const,
%     H = 2(Su' Qbar Su + Rbar),   F = 2 Su' Qbar Sx
%
%   in the form quadprog expects. Only the linear term moves with the
%   state, so each solve costs one Nh x 2 matrix-vector product plus the QP.
%
%   Note this is the NOMINAL error dynamics: the drift d and the reference
%   velocity v_tgt are held constant over the horizon, which is what makes
%   zeta_{k+1} = A zeta_k + B u_k exact rather than approximate. The drift
%   is re-queried and the problem re-solved every inner step, so the
%   variation of d is handled by the receding horizon, not predicted.
%
%   Inputs
%     A,B     : per-axis discrete model (2x2, 2x1)
%     Q,R     : stage weights (2x2, scalar)
%     M       : terminal weight (2x2)
%     Nh      : horizon
%     u_max   : actuator box |u_k| <= u_max (scalar)
%
%   Output (struct qp)
%     qp.H, qp.F   : Hessian and state-to-linear-term map
%     qp.Sx, qp.Su : prediction matrices, kept for testing
%     qp.lb, qp.ub : box bounds, Nh x 1
%     qp.opts      : quadprog options (trust-region-reflective, so the
%                    warm start x0 is actually used; it accepts
%                    bound-only problems, which is exactly this one)

    % A^k for k = 0..Nh, so Apow{k+1} = A^k
    Apow = cell(Nh + 1, 1);
    Apow{1} = eye(size(A));
    for k = 1:Nh
        Apow{k+1} = Apow{k} * A;
    end

    Sx = zeros(2*Nh, 2);
    Su = zeros(2*Nh, Nh);
    for k = 1:Nh
        Sx(2*k-1:2*k, :) = Apow{k+1};
        for j = 1:k
            Su(2*k-1:2*k, j) = Apow{k-j+1} * B;
        end
    end

    % Stage weight on zeta_1..zeta_{Nh-1}, terminal weight on zeta_Nh.
    Qbar = zeros(2*Nh);
    for k = 1:Nh-1
        Qbar(2*k-1:2*k, 2*k-1:2*k) = Q;
    end
    Qbar(2*Nh-1:2*Nh, 2*Nh-1:2*Nh) = M;
    Rbar = R * eye(Nh);

    H = 2 * (Su' * Qbar * Su + Rbar);
    H = (H + H') / 2;               % kill asymmetry from round-off
    F = 2 * (Su' * Qbar * Sx);

    qp = struct();
    qp.H  = H;
    qp.F  = F;
    qp.Sx = Sx;
    qp.Su = Su;
    qp.lb = -u_max * ones(Nh, 1);
    qp.ub =  u_max * ones(Nh, 1);
    qp.opts = optimoptions('quadprog', ...
                           'Algorithm', 'trust-region-reflective', ...
                           'Display', 'off');
end

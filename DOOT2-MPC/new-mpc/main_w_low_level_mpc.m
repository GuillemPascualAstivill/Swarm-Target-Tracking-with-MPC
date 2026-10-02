clc;close all;
%% DOOT vs independently replanned DOOT + low-level MPC (optional wind).
% Run parameters_w.m first. Each branch solves DOOT from its own positions.
% epsilon_reg stays fixed; no density penalty or new ghost-agent behavior.

use_wind = true;
wind_model = 'uniform';  % 'uniform' or 'taylor_green'; used when use_wind=true
wind_U = 0.01;          % Drift speed scale, in the simulation's velocity units
wind_theta = 0;         % Uniform wind direction (radians)
wind_L = 6;             % Taylor-Green spatial period

if ~exist('params', 'var') || ~isfield(params, 'pos_cont0')
    error('Run parameters_w.m before main_w_low_level_mpc.m.');
end
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'MPC'));

%% Settings
warm_start = false;
warm_start_lambda = false;
show_diagnostics = true;
save_figures = false;
final_plot_limits = [-12, 12, -12, 12];  % Display only: identical across wind runs.

n_inner = 10;
Q_x = 1;
Q_u = 35;
% Velocity-tracking weight. Q_v = 0 is the position-only cost the paper's
% appendix quotes and reproduces the earlier baselines exactly. Nonzero Q_v
% prices the velocity error that e_tr and Delta_t charge for.
%
% Swept on no wind (Q_v = 0, 5, 12, 35; damping 0.707, 0.843, 1.000, 1.407):
% velocity tracking RMS falls monotonically (-9%, -17%, -30%), but final L1
% improves only ~1.2% and the cost is steep -- effort +24% and peak
% acceleration +40% at Q_v = 12, and at Q_v = 35 the peak exceeds the DOOT
% proxy (ratio 1.10 vs 0.573 at Q_v = 0), forfeiting the gentler-actuation
% result. Density error here is set by e_alg and e_discr, not e_tr, so
% sharpening the low level pushes on the wrong term. Kept at 0; the knob
% stays for reproducing the sweep.
Q_v = 0;
Nh = 10;
gamma = 1;
% Actuator box. [] disables it, which leaves the low level an unconstrained
% LQR. With a value set, low_level_mode picks how it is enforced:
%   'clamp' - saturate the unconstrained gain after the fact (saturated LQR)
%   'qp'    - solve the condensed box-constrained problem over the horizon
%             and apply its first move; this is what makes it an MPC
% The two coincide while the box is slack, so check the saturation figure
% printed after the run before reading anything into a null difference.
u_max = [];
low_level_mode = 'qp';
heading_speed_tol = 1e-9;
hf_cutoff_hz = 2.5;  % Same cutoff as the old comparison at 10 Hz sampling.
% Density metrics are computed after Wendland smoothing at h_agent, so they
% cannot resolve agent scatter below that scale. Recompute them at finer
% bandwidths to see where, if anywhere, the two branches separate. Keep
% h_agent in the list so the sweep reproduces the headline numbers.
bw_sweep = [0.5, 1.0, params.h_agent];
bw_sweep_W1 = true;  % One Sinkhorn solve per bandwidth per branch; off if slow.

% Does the agent measure the local wind and cancel it in its reference?
% true  - open-loop cancellation; the plant then behaves as predicted, so
%         the re-solve has nothing left to correct (original behaviour).
% false - the agent knows only its own state; the drift is seen only as the
%         difference between where it planned to be and where it is, and is
%         corrected by re-solving. Expect a standing offset, since the
%         prediction model still carries no d.
% What the agent's control law knows about the wind.
% 'measured'  - it reads the true local wind (assumes an anemometer).
% 'none'      - it knows nothing; the drift is visible only as the gap
%               between where it planned to be and where it is.
% 'estimated' - it infers the drift from that same gap and cancels the
%               estimate. Needs only self-localization. NOTE: in a
%               noise-free simulation with exact dynamics the estimator
%               recovers d essentially exactly, so with drift_alpha = 1
%               this collapses onto 'measured'. The difference is in what
%               the agent must sense, not in the numbers.
drift_mode = 'estimated';
drift_alpha = 0.2;   % estimator gain in (0,1]; 1 is deadbeat

drift_feedforward = strcmp(drift_mode, 'measured');   % legacy mirror

% How the low-level position reference is anchored each outer step.
% 'reanchor'  - x_ref starts at the agent's measured position, so any
%               displacement left uncorrected is absorbed and forgotten.
%               Paired with drift_mode 'none' a persistent drift then
%               accumulates without bound; paired with 'estimated' the
%               estimate carries the memory instead, which is option 2.
% 'propagate' - x_ref carries forward at the commanded velocity regardless
%               of where the agent actually is, so displacement persists as
%               an error the feedback can work on. Bounded offset instead
%               of an accumulating one, at the risk of the reference
%               running away if the agent cannot keep up.
%
% The paper (Section V.B) describes 'estimated' + 'reanchor'.
ref_mode = 'reanchor';

N = params.N;
T = params.T;
dt_outer = params.dt_cont;
dt_inner = dt_outer / n_inner;
assert(T >= 2 && dt_outer > 0 && n_inner >= 1 && mod(n_inner,1) == 0, ...
    'Need T >= 2, positive dt_cont, and a positive integer n_inner.');
assert(hf_cutoff_hz > 0 && hf_cutoff_hz < 0.5/dt_inner, ...
    'HF cutoff must lie between zero and the sampling Nyquist frequency.');
mpc = low_level_mpc_setup(dt_inner, Q_x, Q_u, Nh, gamma, u_max, Q_v, ...
                          low_level_mode);
mpc.drift_feedforward = drift_feedforward;
mpc.drift_mode = drift_mode;
% Per-agent drift estimate, only used when drift_mode is 'estimated'. It is
% built from the agent's own motion: the plant adds dt_inner*d to position
% on top of what the drift-free model predicts, so the residual divided by
% dt_inner is a one-step measurement of d, smoothed here by drift_alpha.
d_hat = zeros(N, 2);
% Warm start for the QP path, one horizon per agent per axis. In a receding
% horizon the solution barely moves between inner steps, so feeding the
% shifted previous solution back cuts solver iterations sharply.
U_ws = zeros(Nh, 2, N);
qp_active_total = 0;
qp_calls = 0;
if use_wind
    addpath(fullfile(fileparts(script_dir), 'WINDS'));
    switch lower(wind_model)
        case 'uniform'
            wind_fn = @(tt,p) uniform_wind(tt,p,wind_U,wind_theta);
        case 'taylor_green'
            wind_fn = @(tt,p) taylor_green_wind(tt,p,wind_U,wind_L);
        otherwise
            error('dootMpc:windModel', ...
                'wind_model must be ''uniform'' or ''taylor_green''.');
    end
else
    wind_fn = @(tt,p) zeros(size(p));
end

targetParams = params;
target_rot = isfield(params, 'target_rot') && params.target_rot;
if target_rot
    targetParams.rho_target_sm = [];
end

%% Separate swarm states, solver histories, and measurements
% Reset from the saved initial condition, including when this script is rerun.
names = {'DOOT', 'DOOT+MPC'};
state = struct();
state.name = '';
state.pos = params.pos_cont0;
state.vel = zeros(N,2);
% Anchor the low-level position reference is built from. Under 'reanchor'
% it is overwritten with the measured position every outer step; under
% 'propagate' it advances by dt_outer*v_nom and keeps its own count, so a
% displacement the agent failed to correct survives into the next interval.
state.ref_anchor = params.pos_cont0(1:N,:);
state.alpha = zeros(params.N_tot + params.M,1);
state.lambda = zeros(params.N_tot + params.M,1);
state.pos_hist = zeros(params.N_tot,2,T+1);
state.pos_hist(:,:,1) = state.pos;
state.vel_outer = zeros(T,N,2);
state.vel_inner = zeros(T*n_inner,N,2);
state.cmd_inner = zeros(T*n_inner,N,2);
% DOOT has no physical acceleration state; its entries remain undefined.
state.accel_inner = nan(T*n_inner,N,2);
state.alpha_hist = zeros(params.N_tot+params.M,T);
state.effort = 0;
state.peak_accel = 0;
state.lag_integral = 0;
state.displacement = 0;
state.planned_length = 0;
state.diagn = new_diagnostics(T, params, warm_start, warm_start_lambda);
% Construct the array from a populated template so both elements share fields.
branches = repmat(state,1,2);
for b = 1:2
    branches(b).name = names{b};
end
clear state;

fprintf('DOOT vs DOOT+MPC: two DOOT solves per outer step.\n');
if use_wind
    switch drift_mode
        case 'measured'
            knows = ['MPC measures the local wind and cancels it in its ' ...
                     'reference'];
        case 'none'
            knows = ['MPC does NOT know the wind; it sees the drift only ' ...
                     'through the state it measures, and corrects by ' ...
                     're-solving'];
        case 'estimated'
            knows = sprintf(['MPC does NOT sense the wind; it estimates ' ...
                             'the drift from its own motion (alpha=%g)'], ...
                            drift_alpha);
    end
    fprintf('Wind: %s, U=%g. Both branches are advected; %s.\n', ...
        wind_model, wind_U, knows);
else
    fprintf('Wind: disabled.\n');
end
if strcmp(ref_mode, 'propagate')
    fprintf(['Reference: propagate -- x_ref advances at the commanded ' ...
             'velocity, so an\n  uncorrected displacement survives into ' ...
             'the next interval.\n']);
else
    fprintf(['Reference: reanchor -- x_ref restarts at the measured ' ...
             'position each outer\n  step, so displacement is absorbed ' ...
             'and forgotten.\n']);
end
fprintf('dt_outer=%g, dt_inner=%g, Nh=%d, Q_x=%g, Q_v=%g, Q_u=%g, epsilon_reg=%g\n', ...
    dt_outer, dt_inner, Nh, Q_x, Q_v, Q_u, params.epsilon_reg);
if mpc.constrained
    fprintf('Low level: %s, u_max=%g, N=%d\n', mpc.mode, mpc.u_max, N);
else
    fprintf('Low level: unconstrained LQR, N=%d\n', N);
end
fprintf('MPC gain K = [%.4f %.4f], damping ratio = %.4f\n', ...
    mpc.K(1), mpc.K(2), mpc.K(2)/(2*sqrt(mpc.K(1))));

% Timing. The low-level block is wrapped per inner step rather than per
% agent: at N=500 that is 5e3 tic/toc pairs instead of 2.5e6, so the
% instrumentation itself stays far below the measurement.
timing = struct('doot', 0, 'lowlevel', 0, 'loop', 0, 'other', 0);
t_loop = tic;

for t = 1:T
    % Both branches see the same target, held over this outer interval.
    if target_rot
        targetParams.rho_target = star_target_w(params, params.target_omega*(t-1));
    end
    for b = 1:2
        t_doot = tic;
        [v_doot, alpha, lambda, solver_diag] = plan_doot( ...
            branches(b).pos, branches(b).alpha, branches(b).lambda, ...
            targetParams, warm_start, warm_start_lambda);
        timing.doot = timing.doot + toc(t_doot);
        fields = fieldnames(solver_diag);
        for j = 1:numel(fields)
            branches(b).diagn.(fields{j})(t) = solver_diag.(fields{j});
        end
        branches(b).alpha = alpha;
        branches(b).lambda = lambda;
        branches(b).alpha_hist(:,t) = alpha;
        p_start = branches(b).pos(1:N,:);
        v_previous = branches(b).vel;

        if b == 1
            % First-order DOOT has instantaneous command changes. This is an
            % OUTER-STEP finite-difference effort proxy, not physical effort.
            a_proxy = (v_doot-v_previous)/dt_outer;
            branches(b).effort = branches(b).effort + sum(a_proxy(:).^2)*dt_outer;
            branches(b).peak_accel = max(branches(b).peak_accel, ...
                max(vecnorm(a_proxy,2,2)));
            for s = 1:n_inner
                idx = (t-1)*n_inner+s;
                branches(b).vel_inner(idx,:,:) = reshape(v_doot,[1,N,2]);
                branches(b).cmd_inner(idx,:,:) = reshape(v_doot,[1,N,2]);
                if use_wind
                    time_now = (t-1)*dt_outer+(s-1)*dt_inner;
                    d_local = wind_fn(time_now,branches(b).pos(1:N,:));
                    % DOOT has no wind compensation. Re-query spatial wind
                    % at the same cadence as MPC, with the command held.
                    branches(b).pos(1:N,:) = branches(b).pos(1:N,:) + ...
                        dt_inner*(v_doot+d_local);
                    branches(b).lag_integral = branches(b).lag_integral + ...
                        sum(d_local(:).^2)*dt_inner;
                end
            end
            branches(b).vel = v_doot;
            if ~use_wind
                branches(b).pos(1:N,:) = p_start + dt_outer*v_doot;
            end
        else
            % Where this interval's reference starts. Re-anchoring forgets
            % any displacement the agent has accumulated; propagating keeps
            % its own count so the feedback still sees it.
            if strcmp(ref_mode, 'propagate')
                ref_start = branches(b).ref_anchor;
            else
                ref_start = p_start;
            end
            for s = 1:n_inner
                idx = (t-1)*n_inner+s;
                tau = (s-1)*dt_inner;
                time_now = (t-1)*dt_outer + tau;
                p_ref = ref_start + tau*v_doot;
                % Both spectra use velocities at the START of each interval.
                branches(b).vel_inner(idx,:,:) = reshape(branches(b).vel,[1,N,2]);
                branches(b).cmd_inner(idx,:,:) = reshape(v_doot,[1,N,2]);
                d_local = wind_fn(time_now,branches(b).pos(1:N,:));
                err_v = branches(b).vel+d_local-v_doot;
                branches(b).lag_integral = branches(b).lag_integral + ...
                    sum(err_v(:).^2)*dt_inner;
                t_ll = tic;
                for m = 1:N
                    p_old = branches(b).pos(m,:);
                    v_old = branches(b).vel(m,:);
                    [u, v_next, ~, info] = low_level_mpc(p_old, v_old, p_ref(m,:), ...
                        v_doot(m,:), time_now, wind_fn, mpc, U_ws(:,:,m), ...
                        d_hat(m,:));
                    if ~isempty(info.U_next)
                        U_ws(:,:,m) = info.U_next;
                    end
                    qp_active_total = qp_active_total + info.n_active;
                    qp_calls = qp_calls + 1;
                    branches(b).vel(m,:) = v_next;
                    % Match setup's B=[dt^2/2;dt]. The shared helper's third
                    % output uses semi-implicit Euler, so integrate here.
                    % Wind is held locally for this applied substep; spatial
                    % variation is not predicted over the MPC horizon.
                    branches(b).pos(m,:) = p_old + dt_inner*(v_old+info.d) + ...
                        0.5*dt_inner^2*u;
                    if strcmp(drift_mode, 'estimated')
                        % Residual against the drift-free prediction. The
                        % agent needs only its own successive positions for
                        % this -- no wind sensor.
                        p_pred = p_old + dt_inner*v_old + 0.5*dt_inner^2*u;
                        d_meas = (branches(b).pos(m,:) - p_pred) / dt_inner;
                        d_hat(m,:) = (1-drift_alpha)*d_hat(m,:) + ...
                                     drift_alpha*d_meas;
                    end
                    branches(b).accel_inner(idx,m,:) = reshape(u,[1,1,2]);
                    branches(b).effort = branches(b).effort + sum(u.^2)*dt_inner;
                    branches(b).peak_accel = max(branches(b).peak_accel,norm(u));
                end
                timing.lowlevel = timing.lowlevel + toc(t_ll);
            end
            % Advance the anchor by the commanded motion, independently of
            % where the agent actually ended up.
            branches(b).ref_anchor = ref_start + dt_outer*v_doot;
        end

        branches(b).displacement = branches(b).displacement + ...
            sum(vecnorm(branches(b).pos(1:N,:)-p_start,2,2));
        branches(b).planned_length = branches(b).planned_length + ...
            sum(vecnorm(v_doot,2,2))*dt_outer;
        branches(b).vel_outer(t,:,:) = reshape(branches(b).vel,[1,N,2]);
        branches(b).pos_hist(:,:,t+1) = branches(b).pos;
        speed = vecnorm(branches(b).vel,2,2);
        branches(b).diagn.mean_speed(t) = mean(speed);
        branches(b).diagn.max_speed(t) = max(speed);

        % Rebuild density AFTER movement, compared with this interval's target.
        rho = agent_density(branches(b).pos(1:N,:),params);
        error_rho = rho-targetParams.rho_target;
        branches(b).diagn.err_L1(t) = sum(abs(error_rho))*params.dx*params.dy;
        branches(b).diagn.err_L2(t) = sqrt(sum(error_rho.^2)*params.dx*params.dy);
        branches(b).rho_final = rho;
    end
    if mod(t,10) == 0 || t == T
        fprintf('Step %d/%d: L1 DOOT %.4f, DOOT+MPC %.4f\n', t,T, ...
            branches(1).diagn.err_L1(t), branches(2).diagn.err_L1(t));
    end
end

%% Timing breakdown
timing.loop  = toc(t_loop);
timing.other = timing.loop - timing.doot - timing.lowlevel;
fprintf('\n--- Timing (N=%d, T=%d, n_inner=%d) ---\n', N, T, n_inner);
fprintf('  DOOT solves     %8.1f s  (%5.1f%%)\n', ...
    timing.doot, 100*timing.doot/timing.loop);
fprintf('  Low level       %8.1f s  (%5.1f%%)\n', ...
    timing.lowlevel, 100*timing.lowlevel/timing.loop);
fprintf('  Other in loop   %8.1f s  (%5.1f%%)\n', ...
    timing.other, 100*timing.other/timing.loop);
fprintf('  Loop total      %8.1f s  (%.1f min)\n', ...
    timing.loop, timing.loop/60);
if qp_calls > 0
    per_call = timing.lowlevel / qp_calls;
    fprintf('  Low level: %d agent-steps, %.3f ms each\n', ...
        qp_calls, 1000*per_call);
    % What the same low level would cost at the headline problem size.
    fprintf('  Projected low level at N=500, T=500, n_inner=%d: %.1f min\n', ...
        n_inner, per_call * 500 * 500 * n_inner / 60);
end
fprintf('----------------------------------------\n');

%% Constraint activity -- read this before comparing clamp against QP
if mpc.constrained && qp_calls > 0
    if strcmp(mpc.mode, 'qp')
        denom = qp_calls * Nh * 2;   % every move over every horizon
        what  = 'horizon moves';
    else
        denom = qp_calls * 2;        % only the applied component
        what  = 'applied components';
    end
    fprintf('Actuator box u_max=%g (%s): %d of %d %s saturated (%.4f%%).\n', ...
        mpc.u_max, mpc.mode, qp_active_total, denom, what, ...
        100*qp_active_total/denom);
    if qp_active_total == 0
        fprintf(['  Box never activated -- clamp and QP are identical ' ...
                 'here. Lower u_max for the comparison to mean anything.\n']);
    end
end

%% Metrics: outer-step temporal changes and inner-step spectra
% Exclude the undefined initial heading and rest-to-motion step from turn/
% velocity-change summaries. Effort and peak include startup for both branches.
% Retain main_w's final W1 convention: compare against the smoothed target.
target_final_sm = smooth_on_grid_wendland(targetParams.rho_target,params,params.h_agent);
for b = 1:2
    V = branches(b).vel_outer;
    branches(b).speed = sqrt(sum(V.^2,3));
    branches(b).speed_change = abs(diff(branches(b).speed,1,1));
    branches(b).velocity_change = sqrt(sum(diff(V,1,1).^2,3));
    branches(b).turn = heading_change(V,heading_speed_tol);
    [branches(b).hf, branches(b).frequency, branches(b).power] = ...
        velocity_spectrum(branches(b).vel_inner,dt_inner,hf_cutoff_hz);
    [branches(b).W1, branches(b).W_info] = computeW2( ...
        branches(b).rho_final,target_final_sm,params.x,params.y,'p',1);
    if ~branches(b).W_info.converged
        warning('dootMpc:sinkhorn', '%s: W1 estimate has marginal error %.3e.', ...
            names{b},branches(b).W_info.marg_err);
    end
end
[hf_cmd, f_cmd, power_cmd] = velocity_spectrum( ...
    branches(2).cmd_inner,dt_inner,hf_cutoff_hz);

metric_names = {'Acceleration effort (DOOT proxy)'; 'Peak acceleration (DOOT proxy)'; ...
    'Mean speed'; 'Mean velocity change / outer step'; ...
    'Mean speed change / outer step'; 'Mean absolute heading change (rad)'; ...
    'HF velocity power fraction'; 'Ground velocity tracking RMS (own command)'; ...
    'Displacement / own planned length'; 'Final L1 (raw target)'; ...
    'Final L2 (raw target)'; 'Final W1 estimate (smoothed target)'};
metric_values = zeros(numel(metric_names),2);
for b = 1:2
    metric_values(:,b) = [branches(b).effort; branches(b).peak_accel; ...
        mean(branches(b).speed(:)); mean(branches(b).velocity_change(:)); ...
        mean(branches(b).speed_change(:)); mean(branches(b).turn(:),'omitnan'); ...
        branches(b).hf; sqrt(branches(b).lag_integral/(N*T*dt_outer)); ...
        finite_ratio(branches(b).displacement,branches(b).planned_length); ...
        branches(b).diagn.err_L1(end); branches(b).diagn.err_L2(end); branches(b).W1];
end
ratios = nan(size(metric_values,1),1);
nonzero = abs(metric_values(:,1)) > eps;
ratios(nonzero) = metric_values(nonzero,2)./metric_values(nonzero,1);
comparison = table(metric_names,metric_values(:,1),metric_values(:,2),ratios, ...
    'VariableNames',{'Metric','DOOT','DOOT_MPC','MPC_over_DOOT'});
disp(comparison);
fprintf(['Effort is sum ||u||^2 dt, not physical power. DOOT effort/peak use\n' ...
    'outer-step velocity differences; MPC uses applied acceleration.\n' ...
    'Speed/change/heading/spectrum use velocity excluding wind drift.\n' ...
    'Tracking RMS compares ground velocity (v+d) with each own DOOT command.\n' ...
    'Displacement ratio is not a transport-success score.\n' ...
    'MPC own-command HF fraction: %.4f (cutoff %.2f Hz).\n'],hf_cmd,hf_cutoff_hz);

%% Resolution sweep: the same density metrics at finer smoothing bandwidths
% Same conventions as the headline table: L1/L2 against the raw target, W1
% against the target smoothed at the SAME bandwidth. Only the agent-side
% kernel width changes, so any bandwidth at which the branches separate marks
% the scale of the structure the headline numbers are blind to.
n_bw = numel(bw_sweep);
bw_L1 = nan(n_bw,2); bw_L2 = nan(n_bw,2); bw_W1 = nan(n_bw,2);
for k = 1:n_bw
    hk = bw_sweep(k);
    tgt_sm_k = smooth_on_grid_wendland(targetParams.rho_target,params,hk);
    for b = 1:2
        rho_k = agent_density(branches(b).pos(1:N,:),params,hk);
        err_k = rho_k-targetParams.rho_target;
        bw_L1(k,b) = sum(abs(err_k))*params.dx*params.dy;
        bw_L2(k,b) = sqrt(sum(err_k.^2)*params.dx*params.dy);
        if bw_sweep_W1
            [bw_W1(k,b),info_k] = computeW2( ...
                rho_k,tgt_sm_k,params.x,params.y,'p',1);
            if ~info_k.converged
                warning('dootMpc:sinkhornSweep', ...
                    '%s at h=%.2f: W1 estimate has marginal error %.3e.', ...
                    names{b},hk,info_k.marg_err);
            end
        end
    end
end
rows_per_bw = 2+double(bw_sweep_W1);
sweep_names = cell(n_bw*rows_per_bw,1);
sweep_values = zeros(n_bw*rows_per_bw,2);
for k = 1:n_bw
    base = (k-1)*rows_per_bw;
    sweep_names{base+1} = sprintf('L1 (h=%.2f)',bw_sweep(k));
    sweep_names{base+2} = sprintf('L2 (h=%.2f)',bw_sweep(k));
    sweep_values(base+1,:) = bw_L1(k,:);
    sweep_values(base+2,:) = bw_L2(k,:);
    if bw_sweep_W1
        sweep_names{base+3} = sprintf('W1 (h=%.2f)',bw_sweep(k));
        sweep_values(base+3,:) = bw_W1(k,:);
    end
end
sweep_ratios = nan(size(sweep_values,1),1);
nonzero_sweep = abs(sweep_values(:,1)) > eps;
sweep_ratios(nonzero_sweep) = ...
    sweep_values(nonzero_sweep,2)./sweep_values(nonzero_sweep,1);
resolution = table(sweep_names,sweep_values(:,1),sweep_values(:,2),sweep_ratios, ...
    'VariableNames',{'Metric','DOOT','DOOT_MPC','MPC_over_DOOT'});
disp(resolution);

%% Dispersion: unsmoothed agent statistics
% These use agent positions directly, so no kernel bandwidth can hide the
% scatter. The target references are moments of the raw target density.
w_tgt = targetParams.rho_target(:);
w_tgt = w_tgt/sum(w_tgt);
c_tgt = w_tgt'*params.grid_pts;
spread_tgt = sqrt(w_tgt'*sum((params.grid_pts-c_tgt).^2,2));
nn_mean = zeros(1,2); nn_median = zeros(1,2);
cloud_spread = zeros(1,2); centroid_offset = zeros(1,2); off_grid = zeros(1,2);
for b = 1:2
    P = branches(b).pos(1:N,:);
    d2 = (P(:,1)-P(:,1)').^2+(P(:,2)-P(:,2)').^2;
    d2(1:N+1:end) = Inf;  % Self-pairs are not neighbours.
    nn = sqrt(min(d2,[],2));
    nn_mean(b) = mean(nn);
    nn_median(b) = median(nn);
    c = mean(P,1);
    cloud_spread(b) = sqrt(mean(sum((P-c).^2,2)));
    centroid_offset(b) = norm(c-c_tgt);
    off_grid(b) = mean(P(:,1)<params.x(1) | P(:,1)>params.x(end) | ...
                       P(:,2)<params.y(1) | P(:,2)>params.y(end));
end
dispersion_names = {'Mean nearest-neighbour distance'; ...
    'Median nearest-neighbour distance'; 'Cloud RMS radius about own centroid'; ...
    'Centroid offset from target centroid'; 'Fraction of agents off the grid'};
dispersion_values = [nn_mean; nn_median; cloud_spread; centroid_offset; off_grid];
dispersion_ratios = nan(size(dispersion_values,1),1);
nonzero_disp = abs(dispersion_values(:,1)) > eps;
dispersion_ratios(nonzero_disp) = ...
    dispersion_values(nonzero_disp,2)./dispersion_values(nonzero_disp,1);
dispersion = table(dispersion_names,dispersion_values(:,1), ...
    dispersion_values(:,2),dispersion_ratios, ...
    'VariableNames',{'Metric','DOOT','DOOT_MPC','MPC_over_DOOT'});
disp(dispersion);
fprintf(['Target RMS radius %.4f; target centroid (%.4f, %.4f).\n' ...
    'At bandwidths near or below the mean nearest-neighbour distance the KDE\n' ...
    'resolves individual agents, so L1/L2 levels inflate for BOTH branches\n' ...
    'from finite-sample noise; only the ratio column is meaningful there.\n' ...
    'A nonzero off-grid fraction means mass has left the metric domain and\n' ...
    'every grid-based number above understates that branch error.\n'], ...
    spread_tgt,c_tgt(1),c_tgt(2));

%% Comparison plots (no calls to plotresults_w, which overwrites fixed names)
colors = [0 .45 .74; .85 .33 .10];
outer_time = (1:T)'*dt_outer;
change_time = (2:T)'*dt_outer;
fig_motion = figure('Name','DOOT vs DOOT+MPC: motion','Position',[80 80 1400 800]);
for panel = 1:4
    subplot(2,3,panel); hold on; grid on;
    fields = {'speed','velocity_change','speed_change','turn'};
    labels = {'Speed','Velocity change per outer step', ...
        'Speed change per outer step','Absolute heading change (rad)'};
    if panel == 1, tt = outer_time; else, tt = change_time; end
    h = gobjects(2,1);
    for b = 1:2
        h(b) = plot_ensemble(tt,branches(b).(fields{panel}),colors(b,:));
    end
    xlabel('Time (s)'); ylabel(labels{panel});
    title([labels{panel}, ' (mean; dashed 10/90 percentiles)']);
    legend(h,names,'Location','best');
end
subplot(2,3,5); hold on; grid on;
plot(outer_time,branches(1).diagn.err_L1,'Color',colors(1,:),'LineWidth',1.4);
plot(outer_time,branches(2).diagn.err_L1,'Color',colors(2,:),'LineWidth',1.4);
xlabel('Time (s)'); ylabel('L1 density error'); title('Transport vs interval target');
legend(names,'Location','best');
subplot(2,3,6);
loglog(branches(1).frequency(2:end),max(branches(1).power(2:end),realmin), ...
    'Color',colors(1,:),'LineWidth',1.3); hold on;
loglog(branches(2).frequency(2:end),max(branches(2).power(2:end),realmin), ...
    'Color',colors(2,:),'LineWidth',1.3);
loglog(f_cmd(2:end),max(power_cmd(2:end),realmin),'k--','LineWidth',1);
xlabel('Frequency (Hz)'); ylabel('Velocity power per bin'); grid on;
title('Mean per-agent power spectra');
legend({'DOOT','DOOT+MPC','MPC own DOOT command'},'Location','best');

fig_final = figure('Name','DOOT vs DOOT+MPC: final distribution', ...
    'Position',[100 100 1200 480]);
final_layout = tiledlayout(fig_final,1,3,'TileSpacing','compact','Padding','compact');
fields_final = [targetParams.rho_target,branches(1).rho_final,branches(2).rho_final];
% Use the initial target only: identical limits/levels for every wind mode
% with the same parameters. Clip display colors, never densities or metrics.
density_clim = [0, max(max(params.rho_target(:)),eps)];
density_levels = linspace(density_clim(1),density_clim(2),21);
titles = {'Target (last planning interval)','DOOT','DOOT+MPC'};
for panel = 1:3
    nexttile(final_layout,panel);
    density_display = min(fields_final(:,panel),density_clim(2));
    contourf(params.X,params.Y,reshape(density_display, ...
        size(params.X)),density_levels,'LineColor','none'); hold on;
    if panel > 1
        p = branches(panel-1).pos;
        scatter(p(1:N,1),p(1:N,2),4,'w','filled');
    end
    axis equal;
    axis(final_plot_limits);
    set(gca,'Color',[0.65 0.65 0.65]);  % Keep white agents visible outside the grid.
    caxis(density_clim); colormap(gca,parula);
    cb = colorbar;
    ylabel(cb,'Density');
    title(titles{panel});
end
xlabel(final_layout, { ...
    'Colors: grid-normalized density; white dots: agents. Colors saturate at the common upper limit.', ...
    sprintf('Density is evaluated only on [%g, %g] x [%g, %g]; gray areas outside this grid have no density evaluation.', ...
        min(params.X(:)),max(params.X(:)),min(params.Y(:)),max(params.Y(:)))}, ...
    'FontSize',10,'Interpreter','none');
if show_diagnostics
    for b = 1:2
        plot_diagnostics_w(branches(b).diagn,params,names{b});
    end
end
if save_figures
    exportgraphics(fig_motion,fullfile(script_dir,'doot2_mpc_motion.pdf'), ...
        'ContentType','vector');
    exportgraphics(fig_final,fullfile(script_dir,'doot2_mpc_final.pdf'), ...
        'ContentType','vector');
end

%% Local helpers
function diagn = new_diagnostics(T,params,warm_start,warm_start_lambda)
    fields = {'d_alpha','rel_d_alpha','alpha_inf','max_gradphi', ...
        'frac_violate','mean_speed','max_speed','w_norm','err_L1','err_L2'};
    diagn = struct();
    for j = 1:numel(fields), diagn.(fields{j}) = zeros(T,1); end
    diagn.warm_start = warm_start;
    diagn.warm_start_lambda = warm_start_lambda;
    diagn.pd_iters = params.pd_iters;
    diagn.alpha_cap_cold = params.etap*params.pd_iters;
end

function [v,alpha,lambda,d] = plan_doot(pos,alpha_prev,lambda_prev,params,ws,wsl)
    % Same planner algebra as main_w, evaluated independently for each swarm.
    N = params.N; N_tot = params.N_tot;
    kernels = computeKernels_w(pos,params.ant_cords,params);
    integral = computeTargetint_w(pos,params.ant_cords,params);
    w_agent = zeros(N_tot,1);
    w_agent(1:N) = sum(integral.agents_ag(1:N,1:N),2)/N-integral.target_ag(1:N);
    w_agent(N+1:N_tot) = sum(integral.agents_ag(N+1:N_tot,N+1:N_tot),2)+params.ghost_w;
    if ~isempty(integral.agents_ant)
        int_agents_ant = sum(integral.agents_ant,1)'/N;
        w_ant = params.antw*(int_agents_ant/N-integral.target_ant);
    else
        w_ant = [];
    end
    w_cont = params.w_gain*[w_agent;w_ant];
    w_cont(N+1:N_tot) = w_cont(N+1:N_tot)+params.ghost_w;
    c_sq = conformal(pos,params.ant_cords,params);
    if ws && wsl
        [alpha,lambda] = primal_dual(kernels,integral,w_cont,c_sq,params,alpha_prev,lambda_prev);
    elseif ws
        [alpha,lambda] = primal_dual(kernels,integral,w_cont,c_sq,params,alpha_prev);
    else
        [alpha,lambda] = primal_dual(kernels,integral,w_cont,c_sq,params);
    end
    v = -[kernels.D_x(1:N,:)*alpha,kernels.D_y(1:N,:)*alpha];
    gmag = vecnorm(v,2,2);
    d = struct('d_alpha',norm(alpha-alpha_prev), ...
        'rel_d_alpha',norm(alpha-alpha_prev)/max(norm(alpha_prev),eps), ...
        'alpha_inf',max(abs(alpha)),'max_gradphi',max(gmag), ...
        'frac_violate',mean(gmag>1),'w_norm',norm(w_agent(1:N)));
end

function rho = agent_density(pos,params,h)
    % Chunk the grid to avoid allocating another grid-by-agent-by-axis tensor.
    % h defaults to params.h_agent; pass it explicitly for the resolution sweep.
    if nargin < 3 || isempty(h), h = params.h_agent; end
    rho = zeros(size(params.grid_pts,1),1);
    for first = 1:1000:numel(rho)
        idx = first:min(first+999,numel(rho));
        d2 = (params.grid_pts(idx,1)-pos(:,1)').^2 + ...
             (params.grid_pts(idx,2)-pos(:,2)').^2;
        rho(idx) = mean(wendland_kernel(d2,h),2);
    end
    mass = sum(rho)*params.dx*params.dy;
    if ~isfinite(mass) || mass <= 0
        error('dootMpc:density','Swarm density has no finite positive mass on the grid.');
    end
    rho = rho/mass;  % Match main_w's finite-grid normalization.
end

function theta = heading_change(V,tol)
    a = V(1:end-1,:,:); b = V(2:end,:,:);
    theta = abs(atan2(a(:,:,1).*b(:,:,2)-a(:,:,2).*b(:,:,1),sum(a.*b,3)));
    invalid = sqrt(sum(a.^2,3)) < tol | sqrt(sum(b.^2,3)) < tol;
    theta(invalid) = NaN; % Undefined heading must not count as a zero turn.
end

function [hf,f,P] = velocity_spectrum(V,dt,cutoff_hz)
    V = V-mean(V,1);
    L = size(V,1);
    Y = fft(V,[],1);
    last = floor(L/2)+1;
    P = abs(Y(1:last,:,:)).^2/L^2;
    if mod(L,2) == 0
        P(2:end-1,:,:) = 2*P(2:end-1,:,:);
    else
        P(2:end,:,:) = 2*P(2:end,:,:);
    end
    P = mean(sum(P,3),2);
    f = (0:last-1)'/(L*dt);
    hf = finite_ratio(sum(P(f>=cutoff_hz)),sum(P));
end

function h = plot_ensemble(t,A,color)
    h = plot(t,mean(A,2,'omitnan'),'Color',color,'LineWidth',1.5);
    bounds = prctile(A,[10 90],2);
    plot(t,bounds,'--','Color',color,'LineWidth',0.6,'HandleVisibility','off');
end

function r = finite_ratio(a,b)
    if abs(b) <= eps, r = NaN; else, r = a/b; end
end

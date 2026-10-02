%% run_overnight.m  --  run both paper sweeps back to back (~10.5 h serial)
%
% Before starting:
%   1) set quick_test = false in sweep_wind_seeds.m and sweep_regime_map.m
%   2) stop Windows from sleeping (Settings > System > Power > "Never" when
%      plugged in; also set closing the lid to "Do nothing" on a laptop)
% Then type  run_overnight  in the MATLAB command window and leave it.
% Everything printed goes to overnight_log.txt. If the first sweep fails,
% the second still runs. No variables are used here on purpose: both
% sweeps call parameters_w.m, which clears the workspace.

diary('overnight_log.txt');
fprintf('Started %s\n', char(datetime('now')));

try
    sweep_wind_seeds;           % ~1.5 h -> wind_seeds_results.mat
catch ME
    fprintf(2, 'sweep_wind_seeds FAILED: %s\n', ME.message);
end

try
    sweep_regime_map;           % ~9 h -> regime_map_results.mat (checkpointed)
catch ME
    fprintf(2, 'sweep_regime_map FAILED: %s\n', ME.message);
end

fprintf('Finished %s\n', char(datetime('now')));
diary off;

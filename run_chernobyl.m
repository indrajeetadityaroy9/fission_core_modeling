% RUN_CHERNOBYL - Comprehensive Forensic Simulation
%
% PURPOSE:
% Executes the RBMK-1000 accident simulation with full reactivity breakdown.
% Analyzes ALL contributing factors to the Chernobyl disaster:
%   1. Positive Void Coefficient (the design flaw)
%   2. Doppler Effect (negative feedback from fuel heating)
%   3. Xenon Poisoning (fission product absorption)
%   4. Control Rod Tip Effect (positive SCRAM)
%   5. Spatial Instability (Lower vs Upper core divergence)
%
% OUTPUT:
% Generates forensic timeline tables showing all reactivity components.
%
% DEPENDENCIES:
%   rbmk_parameters.m, rbmk_dynamics.m, compute_equilibrium.m, dde15s_new.m

clear; clc;

%% ========================================================================
%  STEP 1: CONFIGURATION & EQUILIBRIUM
%  ========================================================================
fprintf('============================================================\n');
fprintf(' RBMK-1000 COMPREHENSIVE FORENSIC SIMULATION\n');
fprintf('============================================================\n');

% 1. Load Physics
p = rbmk_parameters();

% 2. Set Accident Conditions
%    - Power: 200 MW (Unstable low-power regime)
%    - Flow:  6000 kg/s (Reduced flow, enhances void formation)
%    - SCRAM: t = 30.0s (AZ-5 button pressed)
target_power_mw = 200;
power_fraction = target_power_mw / 3200;
p.m_flow = 6000;
p.t_scram = 30.0;

fprintf('\n[CONFIG] Accident Conditions:\n');
fprintf('   Target Power:  %6.0f MW (%.1f%% of nominal)\n', target_power_mw, power_fraction*100);
fprintf('   Coolant Flow:  %6.0f kg/s (reduced)\n', p.m_flow);
fprintf('   SCRAM Time:    %6.1f s\n', p.t_scram);

% 3. Criticality Search
[y0, p_sim, diag] = compute_equilibrium(power_fraction, p);

fprintf('\n[STATUS] Initial Equilibrium:\n');
fprintf('   Rod Position (c):    %.4f (%.1f%% inserted)\n', diag.c_eq, diag.c_eq*100);
fprintf('   Void Fraction:       %.2f%%\n', diag.void_fraction * 100);
fprintf('   Fuel Temperature:    %.1f K\n', diag.fuel_temp);
fprintf('   Xenon Concentration: %.2e\n', diag.xenon);
fprintf('   Excess Reactivity:   %+.2f beta\n', diag.rho_0 / p.beta);

fprintf('\n[STATUS] Initial Reactivity Balance:\n');
fprintf('   Void:    %+.4f (%+.2f beta)\n', diag.rho_void, diag.rho_void/p.beta);
fprintf('   Doppler: %+.4f (%+.2f beta)\n', diag.rho_dop, diag.rho_dop/p.beta);
fprintf('   Xenon:   %+.4f (%+.2f beta)\n', diag.rho_xen, diag.rho_xen/p.beta);
fprintf('   Rod:     %+.4f (%+.2f beta)\n', diag.rho_rod, diag.rho_rod/p.beta);
fprintf('   ---------------------------------\n');
fprintf('   Total:   %+.2e (should be ~0)\n', diag.rho_total);

if diag.c_eq > 0.1
    warning('Initial rod position > 0.1. The "Tip Effect" may be diminished.');
else
    fprintf('\n   [OK] Rods withdrawn - Tip effect is ARMED.\n');
end

%% ========================================================================
%  STEP 2: STIFF INTEGRATION
%  ========================================================================
fprintf('\n[SOLVER] Starting dde15s_new integration (0 -> 60s)...\n');

tspan = [0 60];
lags = p_sim.tau_flow;
history = y0;

% Strict tolerances for prompt critical kinetics
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'InitialStep', 1e-4, 'MaxStep', 0.1);

tic;
try
    sol = dde15s_new(@(t,y,Z) rbmk_dynamics(t,y,Z,p_sim), lags, history, tspan, options);
catch ME
    fprintf('\n[ERROR] Solver Failed.\n');
    rethrow(ME);
end
elapsed = toc;
fprintf('[SOLVER] Complete: %.2f seconds, %d steps.\n', elapsed, length(sol.x));

%% ========================================================================
%  STEP 3: EXTRACT ALL STATE VARIABLES
%  ========================================================================
t = sol.x;
y = sol.y;
N = length(t);

% --- LOWER REGION (indices 1-8) ---
n_L   = y(1,:);           % Neutron density (normalized)
C_L   = y(2,:);           % Precursor concentration
alpha_L = y(3,:);         % Void fraction
Tf_L  = y(4,:);           % Fuel temperature (K)
Tm_L  = y(5,:);           % Moderator temperature (K)
I_L   = y(6,:);           % Iodine-135
X_L   = y(7,:);           % Xenon-135
c_L   = y(8,:);           % Control rod position

% --- UPPER REGION (indices 9-16) ---
n_U   = y(9,:);           % Neutron density (normalized)
C_U   = y(10,:);          % Precursor concentration
alpha_U = y(11,:);        % Void fraction
Tf_U  = y(12,:);          % Fuel temperature (K)
Tm_U  = y(13,:);          % Moderator temperature (K)
I_U   = y(14,:);          % Iodine-135
X_U   = y(15,:);          % Xenon-135
c_U   = y(16,:);          % Control rod position

% --- POWER (MW) ---
P_L = n_L * p_sim.k_P;
P_U = n_U * p_sim.k_P;
P_Total = P_L + P_U;

%% ========================================================================
%  STEP 4: CALCULATE ALL REACTIVITY COMPONENTS
%  ========================================================================
% All values in absolute units, then converted to beta

% --- VOID REACTIVITY (Positive - THE RBMK FLAW) ---
rho_void_L = p_sim.kappa_V .* alpha_L;
rho_void_U = p_sim.kappa_V .* alpha_U;

% --- DOPPLER REACTIVITY (Negative - Stabilizing) ---
rho_dop_L = p_sim.kappa_D .* (sqrt(Tf_L) - sqrt(p_sim.T0));
rho_dop_U = p_sim.kappa_D .* (sqrt(Tf_U) - sqrt(p_sim.T0));

% --- XENON REACTIVITY (Negative - Poison) ---
rho_xen_L = -p_sim.kappa_X .* X_L;
rho_xen_U = -p_sim.kappa_X .* X_U;

% --- ROD REACTIVITY (Tip Effect in Lower, Boron only in Upper) ---
% Lower: graphite tip + boron absorption
rho_rod_L = p_sim.kappa_tip .* c_L .* exp(-10.*c_L) - p_sim.kappa_boron .* c_L;
% Upper: boron absorption only (no graphite tip)
rho_rod_U = -p_sim.kappa_boron .* c_U;

% --- EXCESS REACTIVITY (from fuel enrichment) ---
if isfield(p_sim, 'rho_0')
    rho_excess = p_sim.rho_0;
else
    rho_excess = 0;
end

% --- TOTAL REACTIVITY ---
rho_total_L = rho_void_L + rho_dop_L + rho_xen_L + rho_rod_L + rho_excess;
rho_total_U = rho_void_U + rho_dop_U + rho_xen_U + rho_rod_U + rho_excess;

% Convert to beta units for display
beta = p_sim.beta;

%% ========================================================================
%  STEP 5: DEFINE LOG POINTS
%  ========================================================================
% Key moments in the accident sequence
log_times = [0, 10, 20, 29.0, 29.9, 30.0, 30.2, 30.5, 30.7, 31.0, 31.5, 32.0, 33.0, 35.0, 40.0, 50.0, 60.0];

% Find indices closest to log times
log_idx = zeros(size(log_times));
for k = 1:length(log_times)
    [~, log_idx(k)] = min(abs(t - log_times(k)));
end

% Find peak power index
[max_power, peak_idx] = max(P_Total);

%% ========================================================================
%  TABLE 1: POWER & SPATIAL DISTRIBUTION
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' TABLE 1: POWER & SPATIAL DISTRIBUTION\n');
fprintf('============================================================\n');
fprintf('%-7s | %-8s | %-8s | %-8s | %-8s | %s\n', ...
    'Time', 'P_Lower', 'P_Upper', 'P_Total', 'Ratio', 'Phase');
fprintf('  (s)   |   (MW)   |   (MW)   |   (MW)   | (P_L/P_U)|       \n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(log_idx)
    i = log_idx(k);
    ratio = P_L(i) / max(P_U(i), 1e-6);

    % Determine phase
    if t(i) < p_sim.t_scram
        phase = 'Pre-SCRAM';
    elseif t(i) < p_sim.t_scram + 2
        phase = 'Tip Effect';
    elseif t(i) < p_sim.t_scram + 5
        phase = 'Excursion';
    else
        phase = 'Shutdown';
    end

    fprintf('%7.1f | %8.0f | %8.0f | %8.0f | %8.2f | %s\n', ...
        t(i), P_L(i), P_U(i), P_Total(i), ratio, phase);
end

% Print peak separately
fprintf('------------------------------------------------------------\n');
ratio_peak = P_L(peak_idx) / max(P_U(peak_idx), 1e-6);
fprintf('  PEAK  | %8.0f | %8.0f | %8.0f | %8.2f | @ t=%.2fs\n', ...
    P_L(peak_idx), P_U(peak_idx), P_Total(peak_idx), ratio_peak, t(peak_idx));
fprintf('============================================================\n');

%% ========================================================================
%  TABLE 2: REACTIVITY BREAKDOWN (LOWER CORE - WHERE TIP EFFECT OCCURS)
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' TABLE 2: REACTIVITY BREAKDOWN - LOWER CORE (units: beta)\n');
fprintf('============================================================\n');
fprintf('%-7s | %-7s | %-7s | %-7s | %-7s | %-7s | %s\n', ...
    'Time', 'Void', 'Doppler', 'Xenon', 'Rod', 'TOTAL', 'Dominant');
fprintf('  (s)   | (beta)  | (beta)  | (beta)  | (beta)  | (beta)  |         \n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(log_idx)
    i = log_idx(k);

    rv = rho_void_L(i)/beta;
    rd = rho_dop_L(i)/beta;
    rx = rho_xen_L(i)/beta;
    rr = rho_rod_L(i)/beta;
    rt = rho_total_L(i)/beta;

    % Determine dominant effect
    [~, dom_idx] = max(abs([rv, rd, rx, rr]));
    dom_names = {'Void', 'Doppler', 'Xenon', 'Rod'};
    dominant = dom_names{dom_idx};

    fprintf('%7.1f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %s\n', ...
        t(i), rv, rd, rx, rr, rt, dominant);
end

fprintf('------------------------------------------------------------\n');
rv = rho_void_L(peak_idx)/beta;
rd = rho_dop_L(peak_idx)/beta;
rx = rho_xen_L(peak_idx)/beta;
rr = rho_rod_L(peak_idx)/beta;
rt = rho_total_L(peak_idx)/beta;
fprintf('  PEAK  | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | @ t=%.2fs\n', ...
    rv, rd, rx, rr, rt, t(peak_idx));
fprintf('============================================================\n');

%% ========================================================================
%  TABLE 3: REACTIVITY BREAKDOWN (UPPER CORE - NO TIP EFFECT)
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' TABLE 3: REACTIVITY BREAKDOWN - UPPER CORE (units: beta)\n');
fprintf('============================================================\n');
fprintf('%-7s | %-7s | %-7s | %-7s | %-7s | %-7s | %s\n', ...
    'Time', 'Void', 'Doppler', 'Xenon', 'Rod', 'TOTAL', 'Dominant');
fprintf('  (s)   | (beta)  | (beta)  | (beta)  | (beta)  | (beta)  |         \n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(log_idx)
    i = log_idx(k);

    rv = rho_void_U(i)/beta;
    rd = rho_dop_U(i)/beta;
    rx = rho_xen_U(i)/beta;
    rr = rho_rod_U(i)/beta;
    rt = rho_total_U(i)/beta;

    % Determine dominant effect
    [~, dom_idx] = max(abs([rv, rd, rx, rr]));
    dom_names = {'Void', 'Doppler', 'Xenon', 'Rod'};
    dominant = dom_names{dom_idx};

    fprintf('%7.1f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %s\n', ...
        t(i), rv, rd, rx, rr, rt, dominant);
end

fprintf('------------------------------------------------------------\n');
rv = rho_void_U(peak_idx)/beta;
rd = rho_dop_U(peak_idx)/beta;
rx = rho_xen_U(peak_idx)/beta;
rr = rho_rod_U(peak_idx)/beta;
rt = rho_total_U(peak_idx)/beta;
fprintf('  PEAK  | %+7.2f | %+7.2f | %+7.2f | %+7.2f | %+7.2f | @ t=%.2fs\n', ...
    rv, rd, rx, rr, rt, t(peak_idx));
fprintf('============================================================\n');

%% ========================================================================
%  TABLE 4: STATE VARIABLES (TEMPERATURES, VOID, XENON)
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' TABLE 4: STATE VARIABLES - LOWER CORE\n');
fprintf('============================================================\n');
fprintf('%-7s | %-8s | %-8s | %-8s | %-10s | %-7s\n', ...
    'Time', 'Void', 'T_fuel', 'T_mod', 'Xenon', 'Rod Pos');
fprintf('  (s)   |   (%%)    |   (K)    |   (K)    |           |   c    \n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(log_idx)
    i = log_idx(k);
    fprintf('%7.1f | %8.2f | %8.1f | %8.1f | %10.2e | %7.4f\n', ...
        t(i), alpha_L(i)*100, Tf_L(i), Tm_L(i), X_L(i), c_L(i));
end

fprintf('------------------------------------------------------------\n');
fprintf('  PEAK  | %8.2f | %8.1f | %8.1f | %10.2e | %7.4f\n', ...
    alpha_L(peak_idx)*100, Tf_L(peak_idx), Tm_L(peak_idx), X_L(peak_idx), c_L(peak_idx));
fprintf('============================================================\n');

%% ========================================================================
%  TABLE 5: SPATIAL INSTABILITY ANALYSIS
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' TABLE 5: SPATIAL INSTABILITY (LOWER vs UPPER DIVERGENCE)\n');
fprintf('============================================================\n');
fprintf('%-7s | %-9s | %-9s | %-9s | %-9s | %s\n', ...
    'Time', 'rho_L', 'rho_U', 'Delta_rho', 'n_L/n_U', 'Status');
fprintf('  (s)   |  (beta)   |  (beta)   |  (beta)   |           |       \n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(log_idx)
    i = log_idx(k);

    rt_L = rho_total_L(i)/beta;
    rt_U = rho_total_U(i)/beta;
    delta = rt_L - rt_U;
    ratio = n_L(i) / max(n_U(i), 1e-9);

    % Determine stability status
    if abs(delta) < 0.1
        status = 'Symmetric';
    elseif delta > 0
        status = 'L > U';
    else
        status = 'U > L';
    end

    fprintf('%7.1f | %+9.2f | %+9.2f | %+9.2f | %9.2f | %s\n', ...
        t(i), rt_L, rt_U, delta, ratio, status);
end

fprintf('============================================================\n');

%% ========================================================================
%  STEP 6: COMPREHENSIVE VERIFICATION
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' PHYSICS VERIFICATION\n');
fprintf('============================================================\n');

% Check 1: Pre-SCRAM Stability
p_at_29s = interp1(t, P_Total, 29);
if abs(p_at_29s - target_power_mw) < 50
    fprintf('[PASS] Pre-SCRAM Stability: Reactor held %.0f MW (target: %d MW)\n', p_at_29s, target_power_mw);
else
    fprintf('[WARN] Pre-SCRAM Stability: Reactor drifted to %.0f MW\n', p_at_29s);
end

% Check 2: Positive Void Coefficient
max_void_rho = max(rho_void_L)/beta;
if max_void_rho > 0
    fprintf('[PASS] Positive Void Coeff: Peak void reactivity = %+.2f beta\n', max_void_rho);
else
    fprintf('[FAIL] Positive Void Coeff: Void reactivity is negative!\n');
end

% Check 3: Tip Effect (Must go positive after SCRAM)
post_scram_mask = (t > p_sim.t_scram) & (t < p_sim.t_scram + 5);
if any(post_scram_mask)
    max_rod_rho = max(rho_rod_L(post_scram_mask))/beta;
    if max_rod_rho > 0
        fprintf('[PASS] Tip Effect: Peak positive rod reactivity = %+.2f beta\n', max_rod_rho);
    else
        fprintf('[FAIL] Tip Effect: Rod reactivity stayed negative. Check kappa_tip.\n');
    end
end

% Check 4: Doppler Feedback (Should be negative)
min_dop_rho = min(rho_dop_L)/beta;
if min_dop_rho < 0
    fprintf('[PASS] Doppler Feedback: Negative feedback active (min = %+.2f beta)\n', min_dop_rho);
else
    fprintf('[FAIL] Doppler Feedback: Should be negative!\n');
end

% Check 5: Xenon Poisoning (Should be negative)
min_xen_rho = min(rho_xen_L)/beta;
if min_xen_rho < 0
    fprintf('[PASS] Xenon Poisoning: Poison effect present (min = %+.2f beta)\n', min_xen_rho);
else
    fprintf('[FAIL] Xenon Poisoning: Should be negative!\n');
end

% Check 6: Spatial Asymmetry
max_delta_rho = max(abs(rho_total_L - rho_total_U))/beta;
if max_delta_rho > 0.1
    fprintf('[PASS] Spatial Instability: L/U divergence detected (max delta = %.2f beta)\n', max_delta_rho);
else
    fprintf('[INFO] Spatial Instability: Regions remained nearly symmetric\n');
end

% Check 7: Power Excursion
if max_power > 32000
    fprintf('[PASS] Prompt Critical: Power exceeded 32 GW!\n');
elseif max_power > target_power_mw * 2
    fprintf('[WARN] Power Excursion: %.0f MW (%.1fx initial, not prompt critical)\n', ...
        max_power, max_power/target_power_mw);
else
    fprintf('[WARN] No Excursion: Peak power only %.0f MW\n', max_power);
end

fprintf('============================================================\n');

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' FORENSIC SUMMARY\n');
fprintf('============================================================\n');
fprintf('Initial Power:     %8.0f MW\n', P_Total(1));
fprintf('Peak Power:        %8.0f MW @ t = %.2f s\n', max_power, t(peak_idx));
fprintf('Power Multiplier:  %8.1fx\n', max_power / P_Total(1));
fprintf('Time to Peak:      %8.2f s after SCRAM\n', t(peak_idx) - p_sim.t_scram);
fprintf('\nDominant Mechanism at Peak:\n');

% Find dominant at peak
rv = abs(rho_void_L(peak_idx));
rd = abs(rho_dop_L(peak_idx));
rx = abs(rho_xen_L(peak_idx));
rr = abs(rho_rod_L(peak_idx));
[~, dom] = max([rv, rd, rx, rr]);
dom_names = {'Void (positive)', 'Doppler (negative)', 'Xenon (negative)', 'Rod/Tip (positive)'};
fprintf('   %s\n', dom_names{dom});

fprintf('============================================================\n');
fprintf('Simulation Complete.\n');

function run_chernobyl_logs()
% run_chernobyl_logs - Dual-Mode Forensic Simulation
%
% PURPOSE:
%   Executes the RBMK-1000 simulation twice to demonstrate the contrast:
%   1. NORMAL MODE: High power, standard parameters (Safe).
%   2. ACCIDENT MODE: Low power, degraded parameters (Explosion).
%
% OUTPUTS:
%   Detailed forensic logs for both runs for side-by-side comparison.

    clear; clc;
    
    % Load the Master Parameter Struct
    P = rbmk_parameters();
    
    %% SIMULATION 1: NORMAL OPERATION (The Control Case)
    % 3200 MW (100%), Standard Flow, Rods Balanced
    % We force a SCRAM at 30s to prove that it is normally safe.
    p_norm = P.normal;
    p_norm.t_scram = 30.0; 
    
    fprintf('############################################################\n');
    fprintf(' RUN 1: NORMAL OPERATION (High Power, Safe Configuration)\n');
    fprintf('############################################################\n');
    run_simulation(3200, p_norm, 'normal');
    
    fprintf('\n\n');
    
    %% SIMULATION 2: ACCIDENT SCENARIO (The Event)
    % 200 MW (6%), Reduced Flow, Rods Withdrawn, High Burnup
    p_acc = P.accident;
    p_acc.t_scram = 30.0;
    
    fprintf('############################################################\n');
    fprintf(' RUN 2: ACCIDENT SCENARIO (Low Power, Unstable Configuration)\n');
    fprintf('############################################################\n');
    run_simulation(200, p_acc, 'accident');

end

%% ========================================================================
%  CORE SIMULATION ENGINE
%  ========================================================================
function run_simulation(target_mw, p_sim, mode_str)

    %% 1. CONFIGURATION
    pow_frac = target_mw / 3200;
    t_end = 60.0;

    fprintf('[CONFIG] %s Mode Loaded:\n', upper(mode_str));
    fprintf('   Initial Power: %.0f MW (%.1f%%)\n', target_mw, pow_frac*100);
    fprintf('   Flow Rate:     %.0f kg/s\n', p_sim.m_flow);
    fprintf('   Transp Delay:  %.2f s\n', p_sim.tau_flow);
    fprintf('   Neutron Cpl:   %.1f\n', p_sim.Dn);
    fprintf('   SCRAM Time:    %.1f s\n', p_sim.t_scram);

    %% 2. INITIALIZATION
    fprintf('\n[INIT] Calculating Equilibrium...\n');
    [y0, p_sim, diag] = compute_equilibrium(pow_frac, p_sim, mode_str);

    % Metric Check
    fprintf('[METRICS] Initial State Audit:\n');
    fprintf('   Rod Position:      %.4f\n', diag.c_eq);
    fprintf('   Excess Reactivity: %+.2f beta\n', diag.rho_0 / p_sim.beta);
    fprintf('   Void Feedback:     %+.2f beta\n', diag.rho_void / p_sim.beta);
    fprintf('   Doppler Feedback:  %+.2f beta\n', diag.rho_dop / p_sim.beta);

    net_inherent = diag.rho_void + diag.rho_dop;
    if net_inherent > 0
        fprintf('   STATUS: UNSTABLE (Net Positive Feedback: %+.2e)\n', net_inherent);
    else
        fprintf('   STATUS: STABLE (Net Negative Feedback: %+.2e)\n', net_inherent);
    end

    %% 3. EXECUTION
    fprintf('\n[SIM] Running Stiff DDE Integration (0 -> %.0fs)...\n', t_end);
    opts = ddeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 0.05);

    tic;
    try
        sol = dde15s_new(@(t,y,Z) rbmk_dynamics(t,y,Z,p_sim), p_sim.tau_flow, y0, [0 t_end], opts);
    catch ME
        fprintf('\n[CRITICAL FAILURE] Solver Crashed: %s\n', ME.message);
        return;
    end
    fprintf('[SIM] Complete in %.3f seconds. Steps: %d\n', toc, length(sol.x));

    %% 4. FORENSIC ANALYSIS
    t = sol.x;
    y = sol.y;
    N = length(t);

    % Pre-allocate
    P_Tot = zeros(1,N); P_L = zeros(1,N); P_U = zeros(1,N);
    rho_rod = zeros(1,N); rho_void = zeros(1,N); rho_dop = zeros(1,N); rho_net = zeros(1,N);
    spatial_tilt = zeros(1,N);

    for i = 1:N
        % Extract State
        n_L = y(1,i); n_U = y(9,i);
        alpha_L = y(3,i); Tf_L = y(4,i); c_L = y(8,i);

        % Power
        P_L(i) = n_L * 1600; 
        P_U(i) = n_U * 1600;
        P_Tot(i) = P_L(i) + P_U(i);
        spatial_tilt(i) = P_L(i) / max(P_U(i), 1e-6);

        % Reactivity Reconstruction
        beta = p_sim.beta;

        % Rods (Tip Effect)
        rho_rod(i) = (p_sim.kappa_tip * c_L * exp(-10*c_L) - p_sim.kappa_boron * c_L) / beta;

        % Dynamic Coefficients
        n_safe = max(n_L, 0.01);
        
        % Use parameters from p_sim to determine curve
        kv = p_sim.kappa_V_low - (p_sim.kappa_V_low - p_sim.kappa_V_high)*...
             (1 - exp(-n_safe/p_sim.kappa_V_transition));
         
        kd = p_sim.kappa_D_base * (1 - exp(-n_safe * p_sim.doppler_scale));

        rho_void(i) = (kv * alpha_L) / beta;
        rho_dop(i)  = (kd * (sqrt(Tf_L) - sqrt(p_sim.T0))) / beta;

        % Net (Approximate)
        rho_net(i) = rho_rod(i) + rho_void(i) + rho_dop(i) + ...
                     (diag.rho_0/beta) + (diag.rho_xen/beta);
    end

    %% 5. LOGGING TABLES

    % --- TABLE 1: PRE-SCRAM ---
    fprintf('\n========================================================================\n');
    fprintf(' PHASE 1: PRE-SCRAM BEHAVIOR (0s to 30s)\n');
    fprintf('========================================================================\n');
    fprintf('%-8s | %-10s | %-12s | %-12s | %-10s\n', 'Time', 'Power(MW)', 'Void(%)', 'T_Fuel(K)', 'Net Rho($)');
    fprintf('------------------------------------------------------------------------\n');

    check_points = [0, 10, 20, 25, 29, 30];
    for cp = check_points
        idx = find(t >= cp, 1);
        if isempty(idx), continue; end
        fprintf('%6.1fs | %10.1f | %11.2f%% | %11.1f | %+9.3f\n', ...
            t(idx), P_Tot(idx), y(3,idx)*100, y(4,idx), rho_net(idx));
    end

    % --- TABLE 2: SCRAM RESPONSE ---
    fprintf('\n========================================================================\n');
    fprintf(' PHASE 2: SCRAM RESPONSE (AZ-5 Pressed at 30.0s)\n');
    fprintf('========================================================================\n');
    fprintf('%-7s | %-7s | %-10s | %-10s | %-10s | %-8s\n', ...
        'Time', 'Rod(c)', 'Rho_Rod($)', 'Rho_Void($)', 'Power(MW)', 'Tilt(L/U)');
    fprintf('------------------------------------------------------------------------\n');

    scram_window = find(t >= 29.5 & t <= 35.0);
    last_print = -1;
    for ki = 1:length(scram_window)
        k = scram_window(ki);
        if t(k) - last_print > 0.5 || abs(P_Tot(k) - P_Tot(max(1,k-1)))/P_Tot(k) > 0.1
            fprintf('%6.2fs | %7.4f | %+10.3f | %+10.3f | %10.0f | %8.2f\n', ...
                t(k), y(8,k), rho_rod(k), rho_void(k), P_Tot(k), spatial_tilt(k));
            last_print = t(k);
        end
    end

    % --- TABLE 3: OUTCOME ---
    [max_p, idx_p] = max(P_Tot);
    final_p = P_Tot(end);
    
    fprintf('\n========================================================================\n');
    fprintf(' PHASE 3: FINAL OUTCOME\n');
    fprintf('========================================================================\n');
    
    if max_p > target_mw * 2
        % EXPLOSION DETECTED
        fprintf('RESULT: EXPLOSION / POWER EXCURSION\n');
        fprintf('Peak Power:          %.0f MW (%.1fx Nominal)\n', max_p, max_p/3200);
        fprintf('Time of Peak:        %.3f s (%.3fs after SCRAM)\n', t(idx_p), t(idx_p)-p_sim.t_scram);
        
        fprintf('------------------------------------------------------------------------\n');
        fprintf('REACTIVITY AT PEAK ($):\n');
        fprintf('  Rod (Tip):         %+.3f $\n', rho_rod(idx_p));
        fprintf('  Void:              %+.3f $\n', rho_void(idx_p));
        fprintf('  Doppler:           %+.3f $\n', rho_dop(idx_p));
        fprintf('------------------------------------------------------------------------\n');
        
        % Automated Checks
        if rho_rod(idx_p) > 0.5, fprintf('[PASS] Positive Scram Confirmed.\n'); end
        if spatial_tilt(idx_p) > 1.5, fprintf('[PASS] Spatial Decoupling Confirmed.\n'); end
        
    else
        % SAFE SHUTDOWN DETECTED
        fprintf('RESULT: SAFE SHUTDOWN\n');
        fprintf('Power at 60s:        %.1f MW\n', final_p);
        fprintf('Peak Rod Rho:        %+.3f $ (Did not trigger excursion)\n', max(rho_rod(scram_window)));
        
        if final_p < target_mw
            fprintf('[PASS] Reactor shut down successfully.\n');
        else
            fprintf('[FAIL] Reactor failed to shut down.\n');
        end
    end
end

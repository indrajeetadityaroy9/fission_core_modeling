function main()
    % RBMK_BIFURCATION_ANALYSIS Hopf bifurcation study of Chernobyl disaster mechanism
    %
    % Performs bifurcation analysis demonstrating how Hopf bifurcation
    % (thermal-hydraulic oscillations from positive void coefficient + transport
    % delay) creates dangerous low-power instability in the RBMK reactor.
    %
    % Analysis Pipeline:
    %   1. Compute steady-state branch across power range
    %   2. Eigenvalue-based stability analysis via Chebyshev spectral method
    %   3. Detect Hopf bifurcation (stability margin crossing zero)
    %
    % Outputs:
    %   - rbmk_bifurcation_results.mat (all data for post-processing)
    %
    % Total runtime: ~5 minutes on typical hardware
    %
    % Features:
    %   - Pure eigenvalue-based stability analysis (no time-domain simulation)
    %   - Chebyshev pseudospectral method for DDE eigenvalues
    %   - Detailed convergence diagnostics for steady-state solver
    %
    % Note: Visualization functions removed. Use saved .mat file with
    %       external plotting tools if figures are needed.

    fprintf('================================================================\n');
    fprintf('  RBMK BIFURCATION ANALYSIS\n');
    fprintf('  Demonstrating the Chernobyl Disaster Mechanism\n');
    fprintf('================================================================\n\n');

    fprintf('This analysis demonstrates:\n');
    fprintf('  HOPF BIFURCATION - Positive void coefficient + transport delay\n');
    fprintf('                     → Stable fixed point becomes limit cycle at low power\n');
    fprintf('                     → Identifies critical power threshold\n');
    fprintf('                     → Explains why low-power operation was dangerous\n\n');
    fprintf('================================================================\n\n');

    p = rbmk_parameters();

    fprintf('Key parameters (research-calibrated):\n');
    fprintf('  beta    = %.5f (delayed neutron fraction)\n', p.beta);
    fprintf('  kappa_V = %+.3f (void coefficient, POSITIVE)\n', p.kappa_V);
    fprintf('  kappa_D = %+.3f (Doppler coefficient, NEGATIVE)\n', p.kappa_D0);
    fprintf('  tau_flow = %.1f s (transport delay)\n', p.tau_flow);
    fprintf('  rho_tip = %.5f (+1 beta, graphite tip effect)\n\n', p.rho_tip);
    fprintf('================================================================\n\n');

    % Validate dependencies before starting long analysis
    check_dependencies();

    % Start comprehensive logging
    log_filename = sprintf('rbmk_analysis_%s.log', datestr(now, 'yyyymmdd_HHMMSS'));
    diary(log_filename);
    diary on;
    fprintf('================================================================\n');
    fprintf('LOGGING SYSTEM ACTIVATED\n');
    fprintf('Session started: %s\n', datestr(now));
    fprintf('Log file: %s\n', log_filename);
    fprintf('================================================================\n\n');

    total_start_time = tic;

    %% Hopf Bifurcation Analysis
    fprintf('HOPF BIFURCATION ANALYSIS\n');
    fprintf('Sweeping power to identify oscillatory instability threshold...\n\n');

    analysis_time = tic;
    [power_levels, eigenvalues_DDE, stability_margins] = hopf_analysis(p);
    fprintf('\nAnalysis completed in %.1f minutes.\n', toc(analysis_time)/60);

    % Detect Hopf bifurcation from eigenvalue stability margins
    valid_margins = ~isnan(stability_margins);
    if any(valid_margins)
        margins_valid = stability_margins(valid_margins);
        powers_valid = power_levels(valid_margins);

        % Find zero crossing (stable → unstable as power decreases)
        sign_changes = find(diff(sign(margins_valid)) ~= 0);

        if ~isempty(sign_changes)
            idx = sign_changes(1);
            % Linear interpolation to find exact crossing
            P1 = powers_valid(idx);
            P2 = powers_valid(idx+1);
            m1 = margins_valid(idx);
            m2 = margins_valid(idx+1);
            P_hopf = P1 - m1 * (P2 - P1) / (m2 - m1);

            fprintf('\n=> HOPF BIFURCATION at P ≈ %.0f MW (%.1f%% nominal)\n', ...
                P_hopf, 100*P_hopf/p.k_P);

            % Get oscillation frequency at bifurcation
            [~, omega_boundary, ~] = find_hopf_point(power_levels, eigenvalues_DDE, 'Verbose', false);
            if ~isnan(omega_boundary) && omega_boundary > 0.01
                fprintf('   Oscillation frequency: %.4f rad/s (%.3f Hz)\n', omega_boundary, omega_boundary/(2*pi));
                fprintf('   Oscillation period: %.1f s\n', 2*pi/omega_boundary);
            end
        else
            P_hopf = NaN;
            omega_boundary = NaN;
            fprintf('\nWARNING: Could not identify clear Hopf bifurcation point.\n');
        end
    else
        P_hopf = NaN;
        omega_boundary = NaN;
        fprintf('\nWARNING: No valid eigenvalue data available.\n');
    end

    fprintf('================================================================\n\n');

    %% Save Results
    fprintf('SAVING RESULTS\n\n');

    save('rbmk_bifurcation_results.mat', ...
        'power_levels', 'eigenvalues_DDE', 'stability_margins', ...
        'P_hopf', 'p');

    fprintf('Results saved to: rbmk_bifurcation_results.mat\n\n');

    fprintf('================================================================\n\n');

    %% Summary Report
    fprintf('FINAL SUMMARY REPORT\n');
    fprintf('================================================================\n\n');

    total_time = toc(total_start_time);

    fprintf('HOPF BIFURCATION ANALYSIS RESULTS:\n\n');

    if ~isnan(P_hopf)
        fprintf('Critical Power: P_Hopf = %.0f MW (%.1f%% nominal)\n', ...
            P_hopf, 100*P_hopf/p.k_P);
        fprintf('Below threshold: Limit cycle oscillations (UNSTABLE)\n');
        fprintf('Above threshold: Stable fixed point (STABLE)\n');
        fprintf('Mechanism: Positive void feedback + transport delay\n\n');
    else
        fprintf('Could not identify clear bifurcation point\n\n');
    end

    fprintf('CHERNOBYL MECHANISM (Theory):\n');
    fprintf('  - Reactor operated at 200 MW (in Hopf unstable region)\n');
    fprintf('  - System was oscillating due to void-neutron coupling\n');
    fprintf('  - SCRAM triggered, adding +1β reactivity (graphite tips)\n');
    fprintf('  - Oscillation peak + SCRAM exceeded prompt critical (ρ > β)\n');
    fprintf('  - Result: Exponential power divergence → explosion\n\n');

    fprintf('PHYSICAL INTERPRETATION:\n\n');

    fprintf('Why low power is dangerous:\n');
    if ~isnan(P_hopf)
        fprintf('  - Below %.0f MW: System oscillates (Hopf instability)\n', P_hopf);
    end
    fprintf('  - Oscillations create periodic void spikes\n');
    fprintf('  - SCRAM during void spike → rho > beta\n');
    fprintf('  - Weak negative feedbacks (low Doppler, low xenon)\n\n');

    fprintf('Why high power is safer:\n');
    if ~isnan(P_hopf)
        fprintf('  - Above %.0f MW: Stable (no oscillations)\n', P_hopf);
    end
    fprintf('  - Strong Doppler feedback\n');
    fprintf('  - Saturated boiling curve (low void sensitivity)\n');
    fprintf('  - Even if rho > beta briefly, feedbacks contain excursion\n\n');

    fprintf('RBMK Design Flaw:\n');
    fprintf('  - Reactor routinely operated at low power (200 MW)\n');
    fprintf('  - This placed it in the Hopf unstable region\n');
    fprintf('  - Operators unaware that low power meant oscillatory instability\n\n');

    fprintf('================================================================\n\n');

    fprintf('ANALYSIS COMPLETE\n');
    fprintf('Total runtime: %.1f minutes\n\n', total_time/60);

    fprintf('Generated files:\n');
    fprintf('  1. rbmk_bifurcation_results.mat (all analysis data)\n\n');

    fprintf('KEY CONCLUSION:\n');
    fprintf('The Hopf bifurcation analysis identifies the critical power threshold\n');
    fprintf('below which the RBMK reactor becomes oscillatory unstable.\n');
    fprintf('This explains why low-power operation (200 MW) was inherently dangerous.\n');
    fprintf('================================================================\n\n');

    % Stop logging
    fprintf('================================================================\n');
    fprintf('LOGGING SYSTEM SHUTDOWN\n');
    fprintf('Session ended: %s\n', datestr(now));
    fprintf('Complete log saved to: %s\n', log_filename);
    fprintf('================================================================\n');
    diary off;
end

function check_dependencies()
    % Verify all required functions exist before starting analysis

    required_functions = {
        'hopf_analysis', ...
        'compute_equilibrium_branch', ...
        'analyze_branch_stability', ...
        'find_hopf_point', ...
        'chebyshev_eigenvalues'
    };

    missing = {};
    for i = 1:length(required_functions)
        func_name = required_functions{i};
        if ~exist(func_name, 'file')
            missing{end+1} = func_name;
        end
    end

    if ~isempty(missing)
        fprintf('ERROR: Missing required functions:\n');
        for i = 1:length(missing)
            fprintf('  - %s.m\n', missing{i});
        end
        error('main:MissingDependencies', ...
            'Cannot proceed without required functions. See CLAUDE.md for codebase structure.');
    end

    fprintf('All dependencies verified\n\n');
end

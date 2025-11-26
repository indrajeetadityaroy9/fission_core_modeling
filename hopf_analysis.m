function [power_levels, eigenvalues_DDE, stability_margins_DDE] = hopf_analysis(p, varargin)
    % SWEEP_POWER_HOPF - Identifies Hopf bifurcation point in RBMK reactor
    %
    % Uses a 2-step pipeline to analyze the Hopf bifurcation:
    %   1. Compute steady-state branch across power range
    %   2. Analyze linear stability (eigenvalues) at each equilibrium
    %
    % USAGE:
    %   [powers, eigs, margins] = hopf_analysis(p)
    %   [powers, eigs, margins] = hopf_analysis(p, 'power_range', [100 2000])
    %   [powers, eigs, margins] = hopf_analysis(p, 'n_points', 50)
    %
    % INPUTS:
    %   p - Parameter structure from rbmk_parameters
    %
    % OPTIONAL PARAMETERS (Name-Value pairs):
    %   'power_range'  - [P_min, P_max] in MW (default: [100, 2000])
    %   'n_points'     - Number of points along branch (default: 42)
    %   'verbose'      - Display detailed output (default: true)
    %
    % OUTPUTS:
    %   power_levels         - Power levels tested (MW) [n_points × 1]
    %   eigenvalues_DDE      - Cell array of DDE eigenvalues at each power level
    %   stability_margins_DDE- Stability margin (max real eigenvalue) at each power
    %
    % EXAMPLE:
    %   p = rbmk_parameters();
    %   [powers, eigs, margins] = hopf_analysis(p);
    %
    %   % Detect stability boundary
    %   [P_boundary, omega, info] = find_hopf_point(powers, eigs);
    %   fprintf('Hopf bifurcation at %.0f MW\n', P_boundary);
    %
    % See also: compute_equilibrium_branch, analyze_branch_stability, find_hopf_point

    %% Parse inputs
    parser = inputParser;
    addRequired(parser, 'p', @isstruct);
    addParameter(parser, 'power_range', [100, 2000], @(x) isnumeric(x) && length(x) == 2);
    addParameter(parser, 'n_points', 42, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'verbose', true, @islogical);
    parse(parser, p, varargin{:});

    opts = parser.Results;

    %% Display header
    if opts.verbose
        fprintf('\n');
        fprintf('========================================\n');
        fprintf('   HOPF BIFURCATION ANALYSIS PIPELINE\n');
        fprintf('========================================\n\n');
    end

    %% STEP 1: Compute steady-state branch
    if opts.verbose
        fprintf('STEP 1: Computing steady-state branch...\n');
    end

    [power_levels, steady_states, convergence_info] = compute_equilibrium_branch(p, ...
        'power_range', opts.power_range, ...
        'n_points', opts.n_points, ...
        'verbose', opts.verbose);

    %% STEP 2: Analyze stability along branch
    if opts.verbose
        fprintf('\nSTEP 2: Analyzing linear stability...\n');
    end

    [eigenvalues_DDE, stability_margins_DDE, stability_info] = analyze_branch_stability(...
        steady_states, p, ...
        'N', 12, ...
        'verbose', opts.verbose);

    %% Save detailed iteration data to file
    if opts.verbose
        detail_filename = sprintf('hopf_sweep_details_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));
        fid = fopen(detail_filename, 'w');

        % Write header
        fprintf(fid, 'Power(MW)\tMargin_DDE\tDominant_Real\tDominant_Imag\tFreq(Hz)\tPeriod(s)\n');

        % Write data for each power level
        for i = 1:length(power_levels)
            % Extract eigenvalue information
            if ~isempty(eigenvalues_DDE{i}) && ~isnan(stability_margins_DDE(i))
                eig_dom = eigenvalues_DDE{i}(1);
                dom_real = real(eig_dom);
                dom_imag = imag(eig_dom);

                if abs(dom_imag) > 1e-6
                    freq = abs(dom_imag) / (2*pi);
                    period = 2*pi / abs(dom_imag);
                else
                    freq = 0;
                    period = Inf;
                end
            else
                dom_real = NaN;
                dom_imag = NaN;
                freq = NaN;
                period = NaN;
            end

            % Write row
            fprintf(fid, '%.1f\t%+.6f\t%+.6f\t%+.6f\t%.4f\t%.2f\n', ...
                power_levels(i), stability_margins_DDE(i), dom_real, dom_imag, freq, period);
        end

        fclose(fid);
        fprintf('\n📊 Detailed iteration data saved to: %s\n', detail_filename);
    end

    %% STEP 3: Detect and report Hopf bifurcation
    if opts.verbose
        fprintf('\n========================================\n');
        fprintf('   BIFURCATION ANALYSIS RESULTS\n');
        fprintf('========================================\n\n');
    end

    % Eigenvalue-based Hopf detection
    valid_margins = ~isnan(stability_margins_DDE);
    if any(valid_margins)
        margins_valid = stability_margins_DDE(valid_margins);
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

            hopf_power = P1 - m1 * (P2 - P1) / (m2 - m1);

            if opts.verbose
                fprintf('HOPF BIFURCATION DETECTED:\n');
                fprintf('  Critical power: %.0f MW (%.1f%% nominal)\n', hopf_power, 100*hopf_power/p.k_P);
                fprintf('  Crossing interval: [%.0f, %.0f] MW\n', P2, P1);
            end

            % Get oscillation frequency at bifurcation
            [P_boundary, omega_boundary, boundary_info] = find_hopf_point(...
                power_levels, eigenvalues_DDE, 'Verbose', false);

            if opts.verbose && ~isnan(omega_boundary)
                if omega_boundary < 0.01
                    fprintf('  Instability type: Divergent (exponential growth)\n');
                else
                    fprintf('  Instability type: Oscillatory (Hopf)\n');
                    fprintf('  Oscillation frequency: %.4f rad/s (%.3f Hz)\n', omega_boundary, omega_boundary/(2*pi));
                    fprintf('  Oscillation period: %.1f s\n', 2*pi/omega_boundary);
                end
            end
        else
            if opts.verbose
                fprintf('⚠ No Hopf bifurcation detected in tested range\n');
                fprintf('  All stability margins have same sign\n');
            end
        end
    else
        if opts.verbose
            fprintf('⚠ No valid eigenvalue data available\n');
        end
    end

    if opts.verbose
        fprintf('\n========================================\n');
        fprintf('   ANALYSIS COMPLETE\n');
        fprintf('========================================\n\n');
    end
end

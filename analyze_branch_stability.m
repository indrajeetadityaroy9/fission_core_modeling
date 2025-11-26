function [eigenvalues, stability_margins, stability_info] = analyze_branch_stability(steady_states, p, varargin)
    % ANALYZE_STABILITY_BRANCH - Perform linear stability analysis along steady-state branch
    %
    % Computes eigenvalues and stability margins for each equilibrium point using
    % Chebyshev pseudospectral method for rigorous DDE stability analysis.
    %
    % USAGE:
    %   [eigs, margins, info] = analyze_branch_stability(steady_states, p)
    %   [eigs, margins, info] = analyze_branch_stability(steady_states, p, 'N', 15)
    %   [eigs, margins, info] = analyze_branch_stability(steady_states, p, 'verbose', true)
    %
    % INPUTS:
    %   steady_states - Cell array of 16-element steady-state vectors {n_points × 1}
    %                   (from compute_equilibrium_branch)
    %   p             - Parameter structure from rbmk_parameters
    %
    % OPTIONAL PARAMETERS (Name-Value pairs):
    %   'N'       - Number of Chebyshev collocation points (default: 12)
    %   'verbose' - Display progress information (default: true)
    %
    % OUTPUTS:
    %   eigenvalues       - Cell array of eigenvalue vectors {n_points × 1}
    %                       Each cell contains complex eigenvalues sorted by Re(λ)
    %   stability_margins - Vector of stability margins [n_points × 1]
    %                       margin = max(Re(λ)): negative = stable, positive = unstable
    %   stability_info    - Structure array with analysis details [n_points × 1]
    %                       .dominant_eigenvalue - Rightmost eigenvalue
    %                       .oscillatory_modes - Complex conjugate pairs
    %                       .computation_time - Time elapsed (s)
    %                       .J0 - Instantaneous Jacobian (16×16)
    %                       .J_tau - Delayed Jacobian (16×16 sparse)
    %
    % EXAMPLE:
    %   % Complete workflow
    %   p = rbmk_parameters();
    %   [powers, states, ~] = compute_equilibrium_branch(p);
    %   [eigs, margins, info] = analyze_branch_stability(states, p);
    %
    %   % Identify stable vs unstable regions
    %   stable_idx = margins < 0;
    %   fprintf('Stable points: %d / %d\n', sum(stable_idx), length(margins));
    %
    %   % Find stability boundary
    %   [P_boundary, omega, ~] = find_hopf_point(powers, eigs);
    %
    % See also: compute_equilibrium_branch, compute_stability, find_hopf_point

    %% Parse inputs
    parser = inputParser;
    addRequired(parser, 'steady_states', @iscell);
    addRequired(parser, 'p', @isstruct);
    addParameter(parser, 'N', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'verbose', true, @islogical);
    parse(parser, steady_states, p, varargin{:});

    opts = parser.Results;
    n_points = length(steady_states);

    %% Preallocate outputs
    eigenvalues = cell(n_points, 1);
    stability_margins = zeros(n_points, 1);
    stability_info(n_points) = struct('dominant_eigenvalue', [], ...
        'oscillatory_modes', [], 'computation_time', NaN, 'J0', [], 'J_tau', []);

    %% Display header
    if opts.verbose
        fprintf('=== STABILITY ANALYSIS OF STEADY-STATE BRANCH ===\n');
        fprintf('Number of equilibria: %d\n', n_points);
        fprintf('Method: Chebyshev Pseudospectral (N=%d)\n', opts.N);
        fprintf('Expected time: ~%.0f s (%.1f ms/point)\n\n', n_points * 0.07, 70);
        fprintf('Progress: ');
    end

    %% Analyze stability at each point
    t_start = tic;
    failed_count = 0;

    for i = 1:n_points
        if opts.verbose && mod(i, max(1, floor(n_points/10))) == 0
            fprintf('%d%% ', round(100*i/n_points));
        end

        y_ss = steady_states{i};

        % Skip if steady state failed to converge
        if isempty(y_ss)
            eigenvalues{i} = [];
            stability_margins(i) = NaN;
            failed_count = failed_count + 1;
            continue;
        end

        try
            % Compute eigenvalues using pseudospectral method
            [eigs, margin, info] = compute_stability(y_ss, p, ...
                'N', opts.N, 'verbose', false);

            % Store results
            eigenvalues{i} = eigs;
            stability_margins(i) = margin;
            stability_info(i).dominant_eigenvalue = info.dominant_eigenvalue;
            stability_info(i).oscillatory_modes = info.oscillatory_modes;
            stability_info(i).computation_time = info.computation_time;
            stability_info(i).J0 = info.J0;
            stability_info(i).J_tau = info.J_tau;

        catch ME
            eigenvalues{i} = [];
            stability_margins(i) = NaN;
            stability_info(i).dominant_eigenvalue = NaN;
            stability_info(i).oscillatory_modes = [];
            stability_info(i).computation_time = NaN;
            failed_count = failed_count + 1;

            if opts.verbose
                fprintf('\n✗ Stability analysis failed at point %d: %s\n', i, ME.message);
            end
        end
    end

    total_time = toc(t_start);

    %% Save detailed stability data to file
    if opts.verbose
        detail_filename = sprintf('stability_branch_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));
        fid = fopen(detail_filename, 'w');

        % Write header
        fprintf(fid, 'Index\tStability_Margin\tDominant_Real\tDominant_Imag\tNum_Oscillatory_Modes\tComputation_Time(s)\n');

        % Write data for each point
        for i = 1:n_points
            if ~isempty(eigenvalues{i}) && ~isnan(stability_margins(i))
                dom_eig = stability_info(i).dominant_eigenvalue;
                n_osc = size(stability_info(i).oscillatory_modes, 1);
                comp_time = stability_info(i).computation_time;

                fprintf(fid, '%d\t%+.6f\t%+.6f\t%+.6f\t%d\t%.4f\n', ...
                    i, stability_margins(i), real(dom_eig), imag(dom_eig), ...
                    n_osc, comp_time);
            else
                % Failed analysis
                fprintf(fid, '%d\tNaN\tNaN\tNaN\t0\tNaN\n', i);
            end
        end

        fclose(fid);
        fprintf('\n📊 Detailed stability data saved to: %s\n', detail_filename);
    end

    %% Summary
    if opts.verbose
        fprintf('\n\n=== STABILITY ANALYSIS SUMMARY ===\n');
        fprintf('Total time: %.1f s (avg %.1f ms/point)\n', ...
            total_time, 1000*total_time/n_points);

        valid_margins = ~isnan(stability_margins);
        n_stable = sum(stability_margins < 0);
        n_unstable = sum(stability_margins > 0);

        fprintf('Valid analyses: %d / %d (%.1f%%)\n', ...
            sum(valid_margins), n_points, 100*sum(valid_margins)/n_points);
        fprintf('Stable equilibria: %d (%.1f%%)\n', n_stable, 100*n_stable/n_points);
        fprintf('Unstable equilibria: %d (%.1f%%)\n', n_unstable, 100*n_unstable/n_points);

        if failed_count > 0
            fprintf('⚠ %d analyses failed\n', failed_count);
        else
            fprintf('✓ All stability analyses successful\n');
        end

        % Detect bifurcations
        if any(valid_margins)
            sign_changes = find(diff(sign(stability_margins(valid_margins))) ~= 0);
            if ~isempty(sign_changes)
                fprintf('\n%d bifurcation point(s) detected (sign changes in stability margin)\n', ...
                    length(sign_changes));
            end
        end

        fprintf('===================================\n\n');
    end
end

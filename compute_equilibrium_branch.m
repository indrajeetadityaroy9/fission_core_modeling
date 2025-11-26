function [parameter_values, steady_states, convergence_info] = compute_equilibrium_branch(p, varargin)
    % COMPUTE_STEADY_STATE_BRANCH - Compute steady states along a parameter branch
    %
    % Finds equilibrium solutions across a range of parameter values (typically power level)
    % to construct the steady-state branch for bifurcation analysis.
    %
    % USAGE:
    %   [params, states, info] = compute_equilibrium_branch(p)
    %   [params, states, info] = compute_equilibrium_branch(p, 'power_range', [100 2000])
    %   [params, states, info] = compute_equilibrium_branch(p, 'n_points', 50)
    %
    % INPUTS:
    %   p - Parameter structure from rbmk_parameters
    %
    % OPTIONAL PARAMETERS (Name-Value pairs):
    %   'power_range'  - [P_min, P_max] in MW (default: [100, 2000])
    %   'n_points'     - Number of points along branch (default: 42)
    %   'verbose'      - Display progress information (default: true)
    %
    % OUTPUTS:
    %   parameter_values  - Power levels tested (MW) [n_points × 1]
    %   steady_states     - Cell array of steady-state vectors {n_points × 1}
    %                       Each cell contains 16-element equilibrium state
    %   convergence_info  - Structure array with convergence diagnostics [n_points × 1]
    %                       .converged - Boolean convergence flag
    %                       .max_derivative - Maximum |dy/dt| at equilibrium
    %                       .power_error - Power constraint error
    %                       .solve_time - Computation time (s)
    %
    % EXAMPLE:
    %   % Compute steady-state branch for Hopf analysis
    %   p = rbmk_parameters();
    %   [powers, states, info] = compute_equilibrium_branch(p, ...
    %       'power_range', [100 2000], 'n_points', 50);
    %
    %   % Check convergence
    %   converged = [info.converged];
    %   fprintf('Successfully converged: %d / %d\n', sum(converged), length(converged));
    %
    % See also: compute_equilibrium, analyze_branch_stability

    %% Parse inputs
    parser = inputParser;
    addRequired(parser, 'p', @isstruct);
    addParameter(parser, 'power_range', [100, 2000], @(x) isnumeric(x) && length(x) == 2);
    addParameter(parser, 'n_points', 42, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'verbose', true, @islogical);
    parse(parser, p, varargin{:});

    opts = parser.Results;

    %% Setup parameter sweep
    P_min = opts.power_range(1);
    P_max = opts.power_range(2);
    n_points = opts.n_points;

    % Convert to fractional power (normalized by k_P)
    power_frac_min = P_min / p.k_P;
    power_frac_max = P_max / p.k_P;
    power_levels_frac = linspace(power_frac_min, power_frac_max, n_points);

    parameter_values = power_levels_frac * p.k_P;  % MW

    % Preallocate outputs
    steady_states = cell(n_points, 1);
    convergence_info(n_points) = struct('converged', false, 'max_derivative', NaN, ...
        'power_error', NaN, 'solve_time', NaN);

    %% Display header
    if opts.verbose
        fprintf('=== STEADY-STATE BRANCH COMPUTATION ===\n');
        fprintf('Power range: %.0f - %.0f MW\n', P_min, P_max);
        fprintf('Number of points: %d\n', n_points);
        fprintf('Step size: %.1f MW\n\n', (P_max - P_min) / (n_points - 1));
        fprintf('Progress: ');
    end

    %% Compute steady states along branch
    t_start = tic;
    failed_count = 0;
    prev_solution = [];  % For continuation

    for i = 1:n_points
        if opts.verbose && mod(i, max(1, floor(n_points/10))) == 0
            fprintf('%d%% ', round(100*i/n_points));
        end

        try
            % Compute steady state at this power level
            % Use continuation from previous solution if available
            if isempty(prev_solution)
                [y_ss, info] = compute_equilibrium(power_levels_frac(i), p);
            else
                [y_ss, info] = compute_equilibrium(power_levels_frac(i), p, ...
                    'InitialGuess', prev_solution, 'Verbose', false);
            end

            % Store for continuation (use if reasonably converged)
            if info.converged || info.max_derivative < 0.5
                prev_solution = y_ss;
            end

            % Store results
            steady_states{i} = y_ss;
            convergence_info(i).converged = info.converged;
            convergence_info(i).max_derivative = info.max_derivative;
            convergence_info(i).power_error = info.power_error;
            convergence_info(i).solve_time = info.solve_time;

            % Warn if convergence is poor
            if ~info.converged || info.max_derivative > 1e-6
                failed_count = failed_count + 1;
                if opts.verbose
                    fprintf('\n⚠ Poor convergence at %.0f MW (max |dy/dt| = %.2e)\n', ...
                        parameter_values(i), info.max_derivative);
                end
            end

        catch ME
            failed_count = failed_count + 1;
            steady_states{i} = [];
            convergence_info(i).converged = false;
            convergence_info(i).max_derivative = Inf;
            convergence_info(i).power_error = Inf;
            convergence_info(i).solve_time = NaN;

            if opts.verbose
                fprintf('\n✗ Failed at %.0f MW: %s\n', parameter_values(i), ME.message);
            end
        end
    end

    total_time = toc(t_start);

    %% Save detailed convergence data to file
    if opts.verbose
        detail_filename = sprintf('steady_state_branch_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));
        fid = fopen(detail_filename, 'w');

        % Write header
        fprintf(fid, 'Power(MW)\tConverged\tMax_Derivative\tPower_Error\tSolve_Time(s)\n');

        % Write data for each power level
        for i = 1:n_points
            fprintf(fid, '%.1f\t%d\t%.6e\t%.6e\t%.4f\n', ...
                parameter_values(i), ...
                convergence_info(i).converged, ...
                convergence_info(i).max_derivative, ...
                convergence_info(i).power_error, ...
                convergence_info(i).solve_time);
        end

        fclose(fid);
        fprintf('\n📊 Detailed convergence data saved to: %s\n', detail_filename);
    end

    %% Summary
    if opts.verbose
        fprintf('\n\n=== BRANCH COMPUTATION SUMMARY ===\n');
        fprintf('Total time: %.1f s (avg %.2f s/point)\n', total_time, total_time/n_points);
        fprintf('Successfully converged: %d / %d (%.1f%%)\n', ...
            n_points - failed_count, n_points, 100*(n_points - failed_count)/n_points);

        if failed_count > 0
            fprintf('⚠ %d points failed to converge\n', failed_count);
        else
            fprintf('✓ All points converged successfully\n');
        end
        fprintf('===================================\n\n');
    end
end

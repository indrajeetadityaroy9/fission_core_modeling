function visualize_hopf_spiral(results)
    % VISUALIZE_HOPF_SPIRAL - 3D spiral bifurcation diagram for RBMK analysis
    %
    % Creates a 3D visualization showing how reactor stability degrades as
    % power decreases, with spiraling trajectories that expand below the
    % Hopf bifurcation point.
    %
    % INPUTS:
    %   results - Structure from rbmk_bifurcation_results.mat containing:
    %             .power_levels, .stability_margins, .eigenvalues_DDE, .P_hopf
    %
    % The visualization shows:
    %   - X-axis: Power (MW) - the control parameter
    %   - Y-axis: Neutron density perturbation
    %   - Z-axis: Void fraction perturbation
    %   - Color: Green (stable), Yellow (transition), Red (unstable)

    %% Extract data
    if nargin < 1
        % Load results if not provided
        if exist('rbmk_bifurcation_results.mat', 'file')
            results = load('rbmk_bifurcation_results.mat');
        else
            error('No results file found. Run main.m first.');
        end
    end

    power_values = results.power_levels;
    stability_margins = results.stability_margins;
    eigenvalues_DDE = results.eigenvalues_DDE;
    p = results.p;

    % Get Hopf point
    hopf_power = results.P_hopf;

    % Compute steady states for visualization
    fprintf('Computing steady states for visualization...\n');
    n_points = length(power_values);
    steady_states = cell(n_points, 1);

    for i = 1:n_points
        power_frac = power_values(i) / p.k_P;
        try
            [y_ss, ~] = compute_equilibrium(power_frac, p, 'Verbose', false);
            steady_states{i} = y_ss;
        catch
            steady_states{i} = [];
        end
    end
    fprintf('Steady states computed.\n');

    %% Define stability regions dynamically from results
    P_hopf = hopf_power;          % Critical threshold from analysis

    % Find transition bounds from stability margins (where sign changes)
    % Look for crossings in stability margin
    sign_changes = find(diff(sign(stability_margins)) ~= 0);
    if ~isempty(sign_changes)
        % Get power values around the sign change
        idx_cross = sign_changes(1);
        P_transition_lower = min(power_values(idx_cross), power_values(idx_cross + 1));
        P_transition_upper = max(power_values(idx_cross), power_values(idx_cross + 1));
    else
        % Fallback: use small band around Hopf point
        P_transition_lower = P_hopf - 50;
        P_transition_upper = P_hopf + 50;
    end

    fprintf('Dynamic stability regions:\n');
    fprintf('  Hopf point: %.0f MW\n', P_hopf);
    fprintf('  Transition zone: %.0f - %.0f MW\n', P_transition_lower, P_transition_upper);

    %% Create figure
    fig = figure('Name', 'RBMK Hopf Bifurcation - 3D Spiral Diagram', ...
                 'Position', [100, 100, 1200, 800], ...
                 'Color', 'w');

    %% Main 3D plot
    ax = axes('Position', [0.1, 0.15, 0.85, 0.75]);
    hold on;
    grid on;

    % Parameters for spiral generation
    n_spirals = 3;           % Number of spiral turns to show
    n_points_per_spiral = 100;

    %% Select 4 key points for legend: 2 in red (unstable) zone, 2 in green (stable) zone
    red_powers = [200, 400];      % Unstable region (below 579 MW)
    green_powers = [800, 1200];   % Stable region (above 579 MW)
    legend_powers = [red_powers, green_powers];

    legend_indices = zeros(1, 4);
    for k = 1:4
        [~, legend_indices(k)] = min(abs(power_values - legend_powers(k)));
    end

    h_legend = gobjects(1, 4);  % Initialize graphics handles for legend
    legend_labels = cell(1, 4);  % Initialize legend labels

    %% Plot spirals at ALL power levels (like target image)
    for i = 1:length(power_values)
        P = power_values(i);

        % Skip if no valid steady state
        if isempty(steady_states{i})
            continue;
        end

        y_ss = steady_states{i};

        % Extract equilibrium values
        n_eq = y_ss(1) + y_ss(9);      % Total neutron density
        alpha_eq = y_ss(3) + y_ss(11); % Total void fraction

        % Get stability margin (real part of dominant eigenvalue)
        if i <= length(stability_margins)
            margin = stability_margins(i);
            if i <= length(eigenvalues_DDE) && ~isempty(eigenvalues_DDE{i})
                eigs_i = eigenvalues_DDE{i};
                [~, idx] = max(real(eigs_i));
                omega = abs(imag(eigs_i(idx)));
            else
                omega = 0.01;
            end
        else
            margin = 0;
            omega = 0.01;
        end

        % Determine color based on actual stability margin
        if margin < -1e-5
            color = [0.2, 0.8, 0.2];  % Green - stable
            spiral_scale = 0.002;
            decay_rate = abs(margin) * 50;
        elseif margin > 1e-5
            color = [0.9, 0.2, 0.2];  % Red - unstable
            spiral_scale = 0.02 * (1 + abs(margin) * 100);
            decay_rate = margin * 30;
        else
            color = [1.0, 0.8, 0.0];  % Yellow - transition
            spiral_scale = 0.01;
            decay_rate = margin * 20;
        end

        % Generate spiral trajectory
        t = linspace(0, n_spirals * 2 * pi, n_points_per_spiral);

        % Spiral amplitude (grows or decays based on stability)
        if margin > 0
            amplitude = spiral_scale * exp(decay_rate * t / max(t));
        else
            amplitude = spiral_scale * exp(-abs(decay_rate) * t / max(t));
        end
        amplitude = min(amplitude, 0.15);

        % Spiral coordinates
        if omega > 1e-6
            freq = omega * 100;
        else
            freq = 1;
        end

        n_spiral = n_eq + amplitude .* cos(freq * t);
        alpha_spiral = alpha_eq + amplitude .* sin(freq * t);
        P_spiral = P * ones(size(t));

        % Plot spiral
        plot3(P_spiral, n_spiral, alpha_spiral, '-', ...
              'Color', color, 'LineWidth', 1.5, ...
              'HandleVisibility', 'off');

        % Check if this is one of the 4 legend points
        legend_idx = find(legend_indices == i, 1);
        if ~isempty(legend_idx)
            % This is a legend point - plot larger marker and store handle
            h_legend(legend_idx) = plot3(P, n_eq, alpha_eq, 'o', ...
                  'MarkerSize', 10, ...
                  'MarkerFaceColor', color, ...
                  'MarkerEdgeColor', 'k', ...
                  'LineWidth', 1.5);

            % Store label
            if margin > 1e-5
                legend_labels{legend_idx} = sprintf('%.0f MW: n=%.2f, \\alpha=%.2f (Unstable)', P, n_eq, alpha_eq);
            else
                legend_labels{legend_idx} = sprintf('%.0f MW: n=%.2f, \\alpha=%.2f (Stable)', P, n_eq, alpha_eq);
            end
        else
            % Regular point - smaller marker
            plot3(P, n_eq, alpha_eq, 'o', ...
                  'MarkerSize', 6, ...
                  'MarkerFaceColor', color, ...
                  'MarkerEdgeColor', 'k', ...
                  'LineWidth', 0.5, ...
                  'HandleVisibility', 'off');
        end
    end

    %% Add Hopf bifurcation plane
    ylims = ylim;
    zlims = zlim;

    % Vertical plane at Hopf point
    [Y_plane, Z_plane] = meshgrid(linspace(ylims(1), ylims(2), 10), ...
                                   linspace(zlims(1), zlims(2), 10));
    X_plane = P_hopf * ones(size(Y_plane));

    surf(X_plane, Y_plane, Z_plane, ...
         'FaceColor', [0.5, 0.5, 0.5], ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    % Hopf point label
    text(P_hopf, ylims(2), zlims(2), ...
         sprintf('  Hopf Bifurcation\n  P = %d MW', round(P_hopf)), ...
         'FontSize', 12, 'FontWeight', 'bold', ...
         'Color', [0.3, 0.3, 0.3], ...
         'VerticalAlignment', 'top');

    %% Add transition zone shading
    % Light yellow band for transition region
    [Y_trans, Z_trans] = meshgrid(linspace(ylims(1), ylims(2), 2), ...
                                   linspace(zlims(1), zlims(2), 2));

    % Lower bound of transition
    X_trans_low = P_transition_lower * ones(size(Y_trans));
    surf(X_trans_low, Y_trans, Z_trans, ...
         'FaceColor', [1, 0.9, 0.5], ...
         'FaceAlpha', 0.1, ...
         'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    %% Labels and formatting
    xlabel('Power (MW)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Neutron Density (n)', 'FontSize', 14, 'FontWeight', 'bold');
    zlabel('Void Fraction (\alpha)', 'FontSize', 14, 'FontWeight', 'bold');

    title({'RBMK-1000 Hopf Bifurcation: 3D Spiral Diagram', ...
           'Stability Transition as Power Decreases'}, ...
          'FontSize', 16, 'FontWeight', 'bold');

    % Set view angle
    view([-35, 25]);

    % Axis limits
    xlim([0, max(power_values) * 1.05]);

    %% Add legend with the 4 key operating points
    legend(h_legend, legend_labels, ...
           'Location', 'northwest', ...
           'FontSize', 10, ...
           'Interpreter', 'tex');

    %% Add power scale markers
    % Mark key power levels
    key_powers = [200, P_hopf, 700, 1000];
    key_labels = {'Chernobyl (200 MW)', 'Hopf (579 MW)', 'Safety Limit (700 MW)', '1000 MW'};

    for k = 1:length(key_powers)
        if key_powers(k) <= max(power_values)
            plot3([key_powers(k), key_powers(k)], ylims, [zlims(1), zlims(1)], ...
                  '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, ...
                  'HandleVisibility', 'off');
        end
    end

    hold off;

    %% Save figure
    saveas(fig, 'hopf_spiral_3d.png');
    saveas(fig, 'hopf_spiral_3d.fig');
    fprintf('Saved: hopf_spiral_3d.png and hopf_spiral_3d.fig\n');

end

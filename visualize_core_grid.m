function fig = visualize_core_grid(varargin)
%VISUALIZE_CORE_GRID Creates a grid cross-section of the RBMK reactor core
%
%   fig = visualize_core_grid()
%   fig = visualize_core_grid('animate', true)
%   fig = visualize_core_grid('animate', true, 'rods_out', true)
%
%   Displays a cross-section view of the reactor core as a grid.
%   Interleaving vertical rectangles of fuel rods and graphite moderators.
%   Light blue squares represent coolant in fuel rod channels.
%
%   Optional Parameters:
%       'animate'  - If true, runs neutron chain reaction animation (default: false)
%       'rods_out' - If true, all control rods are lifted (no absorption),
%                    disabling automatic control (default: false)
%
%   Output:
%       fig - Figure handle

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'animate', false, @islogical);
    addParameter(p, 'rods_out', false, @islogical);  % Lift all control rods (no absorption)
    parse(p, varargin{:});
    animate = p.Results.animate;
    rods_out = p.Results.rods_out;

    % Grid dimensions (reduced for better chain reaction visibility)
    n_rows = 24;
    moderator_width = 2;   % Narrow moderator strips
    fuel_width = 5;        % Fuel rod strips
    n_fuel_strips = 4;     % Number of fuel rod regions
    n_moderator_strips = 5; % Moderator on edges and between fuel
    n_cols = n_fuel_strips * fuel_width + n_moderator_strips * moderator_width;  % 30 columns

    % Colors
    color_coolant = [0.15 0.25 0.7];   % Dark blue - cold water coolant (fuel rod regions)
    color_moderator = [0.85 0.85 0.85]; % Light gray - graphite moderator
    color_edge = [0.5 0.5 0.5];        % Grid edge color
    color_fuel_boundary = [0.8 0.2 0.1]; % Red - fuel rod boundary
    color_uranium = [0.0 0.0 0.6];     % Dark blue - uranium fuel pellets
    color_control_rod = [0.2 0.2 0.2]; % Dark gray/black - boron control rods

    % Create a map: true = fuel rod (coolant), false = moderator
    fuel_map = false(n_rows, n_cols);

    % Build interleaving pattern: M(2) | F(5) | M(2) | F(5) | M(2) | F(5) | M(2) | F(5) | M(2)
    col = 1;
    fuel_strip_positions = [];  % Track fuel strip start positions for boundaries
    for i = 1:(n_fuel_strips + n_moderator_strips)
        if mod(i, 2) == 1
            % Odd = moderator strip
            col = col + moderator_width;
        else
            % Even = fuel strip
            fuel_strip_positions = [fuel_strip_positions; col];
            fuel_map(:, col:col+fuel_width-1) = true;
            col = col + fuel_width;
        end
    end

    n_fuel_rods = 4;  % 4 fuel rod strips

    % Control rod positions (center of ALL moderator strips)
    % Moderator strips at: cols 1-2, 8-9, 15-16, 22-23, 29-30
    % Center x positions (0-indexed plot coords): 0.5, 7.5, 14.5, 21.5, 28.5
    n_control_rods = 5;
    control_rod_cols = [0.5, 7.5, 14.5, 21.5, 28.5];  % Center x positions
    control_rod_width = 0.8;  % Narrower than moderator cell
    if rods_out
        control_rod_inserted = false(n_control_rods, 1);  % All rods lifted (no absorption)
    else
        control_rod_inserted = true(n_control_rods, 1);   % Start with all inserted
    end

    % Coolant state arrays (for fuel channel cells where fuel_map is true)
    coolant_temp = zeros(n_rows, n_cols);      % Temperature (0 = cold, 1 = boiling)
    coolant_density = ones(n_rows, n_cols);    % Density (1 = liquid, 0 = void/steam)
    coolant_void = false(n_rows, n_cols);      % True when evaporated

    % Create figure (sized for smaller grid)
    fig = figure('Name', 'RBMK-1000 Core Cross-Section', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 800, 600]);

    hold on;
    axis equal;
    axis tight;

    % Draw all grid squares and store patch handles for dynamic updates
    cell_patches = gobjects(n_rows, n_cols);
    for row = 1:n_rows
        for col = 1:n_cols
            % Define square vertices
            x = [col-1, col, col, col-1, col-1];
            y = [row-1, row-1, row, row, row-1];

            % Color based on whether inside fuel rod region
            if fuel_map(row, col)
                cell_color = color_coolant;
            else
                cell_color = color_moderator;
            end

            cell_patches(row, col) = patch(x, y, cell_color, ...
                  'EdgeColor', color_edge, ...
                  'LineWidth', 0.3);
        end
    end

    % Draw fuel rod boundaries (around each fuel strip)
    for i = 1:length(fuel_strip_positions)
        x_left = fuel_strip_positions(i) - 1;  % Convert to 0-indexed plot coords
        y_top = 0;

        rectangle('Position', [x_left, y_top, fuel_width, n_rows], ...
                  'EdgeColor', color_fuel_boundary, ...
                  'LineWidth', 2.5, ...
                  'LineStyle', '-');
    end

    % Draw control rods in moderator strips (boron absorbers)
    h_control_rods = gobjects(n_control_rods, 1);  % Graphics handles for dynamic updates
    for i = 1:n_control_rods
        x_left = control_rod_cols(i) - control_rod_width/2;
        h_control_rods(i) = rectangle('Position', [x_left, 0, control_rod_width, n_rows], ...
                                       'FaceColor', color_control_rod, ...
                                       'EdgeColor', 'k', ...
                                       'LineWidth', 1);
        % Hide rod if rods_out mode (all rods lifted)
        if rods_out
            set(h_control_rods(i), 'Visible', 'off');
        end
    end

    % Add uranium circles to all fuel rod cells
    % Grey circles = U-238 (fertile material, can absorb neutrons to become Pu-239)
    % Dark blue circles = U-235 (active fissile material)

    rng(42);  % Fixed seed for reproducibility
    u235_fraction = 0.35;     % 35% of fuel cells contain U-235 (high enrichment for demo)
    pellet_radius = 0.35;     % Circle radius within unit square
    color_u238 = [0.5 0.5 0.5];  % Grey - U-238 (fertile material)

    % Get all fuel cell positions
    [fuel_rows, fuel_cols] = find(fuel_map);
    n_fuel_cells = length(fuel_rows);
    n_u235 = round(u235_fraction * n_fuel_cells);

    % Randomly select cells for U-235
    u235_indices = randperm(n_fuel_cells, n_u235);
    u235_set = false(n_fuel_cells, 1);
    u235_set(u235_indices) = true;

    % Draw circles for all fuel cells
    theta = linspace(0, 2*pi, 30);
    for i = 1:n_fuel_cells
        row = fuel_rows(i);
        col = fuel_cols(i);

        % Center of the cell
        cx = col - 0.5;
        cy = row - 0.5;

        % Draw circle
        rx = cx + pellet_radius * cos(theta);
        ry = cy + pellet_radius * sin(theta);

        if u235_set(i)
            % U-235 (active fissile material) - dark blue
            patch(rx, ry, color_uranium, ...
                  'EdgeColor', color_uranium * 0.7, ...
                  'LineWidth', 0.5);
        else
            % U-238 (fertile material) - grey
            patch(rx, ry, color_u238, ...
                  'EdgeColor', color_u238 * 0.7, ...
                  'LineWidth', 0.5);
        end
    end

    % Extract U-235 cell centers for collision detection (used in animation)
    u235_cx = zeros(n_u235, 1);
    u235_cy = zeros(n_u235, 1);
    u235_idx = 1;
    for i = 1:n_fuel_cells
        if u235_set(i)
            u235_cx(u235_idx) = fuel_cols(i) - 0.5;
            u235_cy(u235_idx) = fuel_rows(i) - 0.5;
            u235_idx = u235_idx + 1;
        end
    end

    % Configure axes
    xlim([0, n_cols]);
    ylim([0, n_rows]);

    % Flip y-axis so row 1 is at top
    set(gca, 'YDir', 'reverse');

    % Labels
    xlabel('Column', 'FontSize', 12);
    ylabel('Row', 'FontSize', 12);
    title(sprintf('RBMK-1000 Reactor Core Cross-Section (%d × %d Grid)', n_rows, n_cols), ...
          'FontSize', 14, 'FontWeight', 'bold');

    % Add tick marks
    set(gca, 'XTick', 0:5:n_cols, 'YTick', 0:6:n_rows);
    set(gca, 'TickDir', 'out');
    set(gca, 'FontSize', 10);

    % Add legend
    h_coolant = patch(NaN, NaN, color_coolant, 'EdgeColor', color_edge);
    h_moderator = patch(NaN, NaN, color_moderator, 'EdgeColor', color_edge);
    h_u235 = patch(NaN, NaN, color_uranium, 'EdgeColor', color_uranium * 0.7);
    h_u238 = patch(NaN, NaN, color_u238, 'EdgeColor', color_u238 * 0.7);
    h_boundary = plot(NaN, NaN, '-', 'Color', color_fuel_boundary, 'LineWidth', 2.5);
    h_control_legend = patch(NaN, NaN, color_control_rod, 'EdgeColor', 'k');

    if animate
        % Include neutrons in legend when animating
        h_neutron_legend = scatter(NaN, NaN, 30, 'k', 'filled');
        legend([h_boundary, h_coolant, h_moderator, h_u235, h_u238, h_control_legend, h_neutron_legend], ...
               {sprintf('Fuel Rod Boundary (%d rods)', n_fuel_rods), ...
                'Light Water Coolant', ...
                'Graphite Moderator', ...
                sprintf('U-235 Fissile (%d)', n_u235), ...
                sprintf('U-238 Fertile (%d)', n_fuel_cells - n_u235), ...
                sprintf('Control Rods (%d, boron)', n_control_rods), ...
                'Neutrons'}, ...
               'Location', 'northeastoutside', ...
               'FontSize', 9);
    else
        legend([h_boundary, h_coolant, h_moderator, h_u235, h_u238, h_control_legend], ...
               {sprintf('Fuel Rod Boundary (%d rods)', n_fuel_rods), ...
                'Light Water Coolant', ...
                'Graphite Moderator', ...
                sprintf('U-235 Fissile (%d)', n_u235), ...
                sprintf('U-238 Fertile (%d)', n_fuel_cells - n_u235), ...
                sprintf('Control Rods (%d, boron)', n_control_rods)}, ...
               'Location', 'northeastoutside', ...
               'FontSize', 10);
    end

    % ========== NEUTRON CHAIN REACTION ANIMATION ==========
    if animate
        fprintf('\nStarting neutron chain reaction animation...\n');

        % Animation parameters
        dt = 0.05;                    % Time step per frame
        pause_time = 0.05;            % Pause between frames (slower = more visible)
        neutron_speed = 2.5;          % Neutron speed (grid units/sec)
        collision_radius = 0.35;      % Same as pellet_radius
        max_neutrons = 150;           % Total neutron count limit for balance (increased)
        n_frames = 800;               % Total animation frames

        % Physics parameters
        neutron_threshold = 80;       % Control rod threshold for automatic control (increased)

        % Coolant temperature parameters
        temp_baseline = 0.0;          % Cold water starting point
        temp_evaporate = 0.92;        % Temperature at which water becomes void (high - requires sustained heating)
        temp_condense = 0.25;         % Temperature below which void condenses back
        heat_per_neutron = 0.03;      % Temperature increase per neutron interaction (gradual heating)
        cooling_rate = 0.002;         % Background cooling per frame
        max_coolant_absorption = 0.08; % Maximum absorption probability for cold water (reduced for demo)

        % Neutron balance tracking
        total_produced = 1;           % Start with 1 initial neutron
        total_absorbed_rods = 0;      % Absorbed by control rods
        total_absorbed_coolant = 0;   % Absorbed by water coolant
        total_fissions = 0;           % Fission events
        total_escaped = 0;            % Escaped boundary

        % Neutron state arrays
        neutron_x = [];
        neutron_y = [];
        neutron_vx = [];
        neutron_vy = [];

        % Create scatter plot for neutrons
        h_neutrons = scatter([], [], 30, 'k', 'filled');

        % Start with a single neutron aimed at a random U-235 cell
        target_idx = randi(n_u235);             % Pick a random U-235 as target
        target_x = u235_cx(target_idx);
        target_y = u235_cy(target_idx);

        % Start just inside the first fuel rod (past the left moderator strip)
        % Use column vectors for all neutron state arrays
        neutron_x = [moderator_width + 0.5];    % Start inside first fuel rod
        neutron_y = [target_y];                 % Align with target U-235 vertically

        % Calculate direction toward target
        dx = target_x - neutron_x;
        dy = target_y - neutron_y;
        dist_to_target = sqrt(dx^2 + dy^2);

        if dist_to_target > 0
            neutron_vx = [neutron_speed * (dx / dist_to_target)];
            neutron_vy = [neutron_speed * (dy / dist_to_target)];
        else
            % Target is in same column, move right
            neutron_vx = [neutron_speed];
            neutron_vy = [0];
        end

        fprintf('Initial neutron starting in fuel rod, aimed at U-235 at (%.1f, %.1f)\n', target_x, target_y);
        if rods_out
            fprintf('Control rods: ALL %d LIFTED (rods_out mode - no absorption)\n', n_control_rods);
        else
            fprintf('Control rods: %d total, threshold: %d neutrons\n', n_control_rods, neutron_threshold);
        end

        % Main animation loop
        for frame = 1:n_frames
            % Step 0: Automatic control logic - adjust control rods based on neutron count
            % (disabled when rods_out mode is active - all rods stay lifted)
            if ~rods_out
                active_neutron_count = length(neutron_x);

                if active_neutron_count < neutron_threshold
                    % RAISE every second control rod (rods 2 and 4) to allow more fission
                    for j = 2:2:n_control_rods
                        if control_rod_inserted(j)
                            control_rod_inserted(j) = false;
                            set(h_control_rods(j), 'Visible', 'off');
                            fprintf('Rod %d RAISED (neutrons: %d < %d)\n', j, active_neutron_count, neutron_threshold);
                        end
                    end
                elseif active_neutron_count > neutron_threshold
                    % INSERT all control rods to suppress chain reaction
                    for j = 1:n_control_rods
                        if ~control_rod_inserted(j)
                            control_rod_inserted(j) = true;
                            set(h_control_rods(j), 'Visible', 'on');
                            fprintf('Rod %d INSERTED (neutrons: %d > %d)\n', j, active_neutron_count, neutron_threshold);
                        end
                    end
                end
            end

            % Step 0.5: Background cooling for all coolant cells
            for row = 1:n_rows
                for col = 1:n_cols
                    if fuel_map(row, col)
                        % Gradual cooling toward baseline
                        coolant_temp(row, col) = max(temp_baseline, ...
                            coolant_temp(row, col) - cooling_rate);

                        % Check for condensation (void → liquid)
                        if coolant_void(row, col) && coolant_temp(row, col) < temp_condense
                            coolant_void(row, col) = false;
                            coolant_density(row, col) = 1.0;
                        end

                        % Update density based on temperature (gradual reduction as heating)
                        if ~coolant_void(row, col)
                            coolant_density(row, col) = max(0, 1 - coolant_temp(row, col));
                        end
                    end
                end
            end

            % Step 1: Move neutrons (smooth sub-step movement)
            neutron_x = neutron_x + neutron_vx * dt;
            neutron_y = neutron_y + neutron_vy * dt;

            % Step 2: Check collisions and handle physics
            new_x = []; new_y = []; new_vx = []; new_vy = [];
            remove = false(size(neutron_x));

            for i = 1:length(neutron_x)
                % Check if out of bounds - remove neutron (escaped)
                if neutron_x(i) < 0 || neutron_x(i) > n_cols || ...
                   neutron_y(i) < 0 || neutron_y(i) > n_rows
                    remove(i) = true;
                    total_escaped = total_escaped + 1;
                    continue;
                end

                % Check collision with control rods (100% absorption if inserted)
                absorbed_by_rod = false;
                for j = 1:n_control_rods
                    if control_rod_inserted(j)
                        rod_x = control_rod_cols(j);
                        rod_half_width = control_rod_width / 2;
                        if neutron_x(i) >= rod_x - rod_half_width && ...
                           neutron_x(i) <= rod_x + rod_half_width
                            remove(i) = true;
                            total_absorbed_rods = total_absorbed_rods + 1;
                            absorbed_by_rod = true;
                            break;
                        end
                    end
                end
                if absorbed_by_rod
                    continue;
                end

                % Get current cell (1-indexed)
                cell_col = max(1, min(n_cols, floor(neutron_x(i)) + 1));
                cell_row = max(1, min(n_rows, floor(neutron_y(i)) + 1));

                % Check if in fuel channel (coolant) - heat water and maybe absorb
                if fuel_map(cell_row, cell_col) && ~coolant_void(cell_row, cell_col)
                    % Neutron is in liquid coolant - heat the water
                    coolant_temp(cell_row, cell_col) = min(1.0, ...
                        coolant_temp(cell_row, cell_col) + heat_per_neutron);

                    % Check for evaporation
                    if coolant_temp(cell_row, cell_col) >= temp_evaporate
                        coolant_void(cell_row, cell_col) = true;
                        coolant_density(cell_row, cell_col) = 0;
                    end

                    % Absorption probability based on current density
                    absorption_prob = max_coolant_absorption * coolant_density(cell_row, cell_col);
                    if rand() < absorption_prob
                        remove(i) = true;
                        total_absorbed_coolant = total_absorbed_coolant + 1;
                        continue;
                    end
                end

                % Check collision with U-235 (any neutron can cause fission)
                dist = sqrt((u235_cx - neutron_x(i)).^2 + (u235_cy - neutron_y(i)).^2);
                hit_idx = find(dist < collision_radius, 1);

                % Only allow fission if under neutron limit (balance system)
                current_count = length(neutron_x) - sum(remove) + length(new_x);
                if ~isempty(hit_idx) && current_count < max_neutrons
                    % Fission: spawn 3 neutrons
                    [nx, ny, nvx, nvy] = spawn_neutrons(u235_cx(hit_idx), u235_cy(hit_idx), neutron_speed, 3);
                    new_x = [new_x; nx];
                    new_y = [new_y; ny];
                    new_vx = [new_vx; nvx];
                    new_vy = [new_vy; nvy];
                    remove(i) = true;  % Remove triggering neutron
                    total_fissions = total_fissions + 1;
                    total_produced = total_produced + 3;
                end
            end

            % Remove absorbed/escaped neutrons and add new ones
            % Ensure column vectors using (:) to avoid indexing issues
            keep_idx = find(~remove);
            neutron_x = [neutron_x(keep_idx); new_x(:)];
            neutron_y = [neutron_y(keep_idx); new_y(:)];
            neutron_vx = [neutron_vx(keep_idx); new_vx(:)];
            neutron_vy = [neutron_vy(keep_idx); new_vy(:)];

            % Step 3: Update coolant cell colors based on temperature
            for row = 1:n_rows
                for col = 1:n_cols
                    if fuel_map(row, col)
                        if coolant_void(row, col)
                            % Void cell - make transparent/white
                            set(cell_patches(row, col), 'FaceColor', [1 1 1], 'FaceAlpha', 0.3);
                        else
                            % Liquid coolant - color based on temperature
                            new_color = temp_to_color(coolant_temp(row, col));
                            set(cell_patches(row, col), 'FaceColor', new_color, 'FaceAlpha', 1.0);
                        end
                    end
                end
            end

            % Step 4: Update neutron graphics
            if ~isempty(neutron_x)
                set(h_neutrons, 'XData', neutron_x, 'YData', neutron_y);
            else
                set(h_neutrons, 'XData', [], 'YData', []);
            end
            drawnow;
            pause(pause_time);

            % End animation if no neutrons left - restart with new neutron
            if isempty(neutron_x)
                fprintf('All neutrons lost. Restarting with new neutron...\n');
                target_idx = randi(n_u235);
                target_x = u235_cx(target_idx);
                target_y = u235_cy(target_idx);

                neutron_x = [moderator_width + 0.5];  % Start inside first fuel rod
                neutron_y = [target_y];

                dx = target_x - neutron_x;
                dy = target_y - neutron_y;
                dist_to_target = sqrt(dx^2 + dy^2);

                neutron_vx = [neutron_speed * (dx / dist_to_target)];
                neutron_vy = [neutron_speed * (dy / dist_to_target)];

                total_produced = total_produced + 1;  % Track restart neutron
            end

            % Progress indicator every 100 frames with balance info
            if mod(frame, 100) == 0
                n_voids = sum(coolant_void(:));
                fprintf('Frame %d/%d | Active: %d | Coolant: %d | Rods: %d | Voids: %d | Fissions: %d\n', ...
                        frame, n_frames, length(neutron_x), total_absorbed_coolant, total_absorbed_rods, n_voids, total_fissions);
            end
        end

        % Final balance summary
        n_voids = sum(coolant_void(:));
        fprintf('\n===== NEUTRON BALANCE SUMMARY =====\n');
        fprintf('Total produced:      %d\n', total_produced);
        fprintf('Absorbed (coolant):  %d\n', total_absorbed_coolant);
        fprintf('Absorbed (rods):     %d\n', total_absorbed_rods);
        fprintf('Total escaped:       %d (left boundary)\n', total_escaped);
        fprintf('Total fissions:      %d\n', total_fissions);
        fprintf('Active remaining:    %d\n', length(neutron_x));
        fprintf('Void cells:          %d / %d (%.1f%%)\n', n_voids, n_fuel_cells, 100*n_voids/n_fuel_cells);
        fprintf('Neutron limit:       %d (max allowed)\n', max_neutrons);
        fprintf('===================================\n');
        fprintf('Animation complete.\n');
    end

    hold off;

    fprintf('Core cross-section created: %d x %d grid\n', n_rows, n_cols);
    fprintf('Fuel rod strips: %d (width: %d), Moderator strips: %d (width: %d)\n', ...
            n_fuel_strips, fuel_width, n_moderator_strips, moderator_width);
end

% ========== HELPER FUNCTIONS ==========

function [x, y, vx, vy] = spawn_neutrons(cx, cy, speed, count)
%SPAWN_NEUTRONS Create neutrons at a position with random directions
%   cx, cy - Center position (fission site)
%   speed  - Neutron speed (grid units per second)
%   count  - Number of neutrons to spawn
%
%   Neutrons are placed slightly outside the collision radius (0.35) to
%   prevent immediate re-collision with the atom that was just hit.
%
%   Returns column vectors of positions and velocities

    spawn_offset = 0.4;  % Slightly larger than collision_radius (0.35)

    angles = 2 * pi * rand(count, 1);

    % Place neutrons at spawn_offset distance from center, in direction of travel
    x = cx + spawn_offset * cos(angles);
    y = cy + spawn_offset * sin(angles);

    % Velocity in the same direction as spawn offset
    vx = speed * cos(angles);
    vy = speed * sin(angles);
end

function color = temp_to_color(temp)
%TEMP_TO_COLOR Maps coolant temperature to blue-red spectrum
%   temp - Temperature value (0 = cold, 1 = boiling)
%
%   Color spectrum (blue to red):
%   Cold (0.0)    = Dark blue [0.15, 0.25, 0.7]
%   Cool (0.33)   = Light blue [0.5, 0.7, 1.0]
%   Warm (0.67)   = Light red [1.0, 0.6, 0.5]
%   Hot (1.0)     = Dark red [0.7, 0.15, 0.15]

    if temp < 0.33
        % Cold: dark blue → light blue
        t = temp / 0.33;
        color = (1-t) * [0.15, 0.25, 0.7] + t * [0.5, 0.7, 1.0];
    elseif temp < 0.67
        % Transition: light blue → light red
        t = (temp - 0.33) / 0.34;
        color = (1-t) * [0.5, 0.7, 1.0] + t * [1.0, 0.6, 0.5];
    else
        % Hot: light red → dark red
        t = min(1, (temp - 0.67) / 0.33);
        color = (1-t) * [1.0, 0.6, 0.5] + t * [0.7, 0.15, 0.15];
    end
end
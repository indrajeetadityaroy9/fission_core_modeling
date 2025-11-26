function fig = visualize_core_grid(varargin)
%VISUALIZE_CORE_GRID Creates a 56x56 grid cross-section of the RBMK reactor core
%
%   fig = visualize_core_grid()
%   fig = visualize_core_grid('animate', true)
%
%   Displays a cross-section view of the reactor core as a 56x56 grid.
%   Interleaving vertical rectangles of fuel rods and graphite moderators.
%   Light blue squares represent coolant in fuel rod channels.
%
%   Optional Parameters:
%       'animate' - If true, runs neutron chain reaction animation (default: false)
%
%   Output:
%       fig - Figure handle

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'animate', false, @islogical);
    parse(p, varargin{:});
    animate = p.Results.animate;

    % Grid dimensions (reduced for better chain reaction visibility)
    n_rows = 24;
    moderator_width = 2;   % Narrow moderator strips
    fuel_width = 5;        % Fuel rod strips
    n_fuel_strips = 4;     % Number of fuel rod regions
    n_moderator_strips = 5; % Moderator on edges and between fuel
    n_cols = n_fuel_strips * fuel_width + n_moderator_strips * moderator_width;  % 30 columns

    % Colors
    color_coolant = [0.7 0.85 1.0];    % Light blue - water coolant (fuel rod regions)
    color_moderator = [0.85 0.85 0.85]; % Light gray - graphite moderator
    color_edge = [0.5 0.5 0.5];        % Grid edge color
    color_fuel_boundary = [0.8 0.2 0.1]; % Red - fuel rod boundary
    color_uranium = [0.0 0.0 0.6];     % Dark blue - uranium fuel pellets

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

    % Create figure (sized for smaller grid)
    fig = figure('Name', 'RBMK-1000 Core Cross-Section', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 800, 600]);

    hold on;
    axis equal;
    axis tight;

    % Draw all grid squares
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

            patch(x, y, cell_color, ...
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

    % Add uranium circles to all fuel rod cells
    % Grey circles = U-238 (fertile material, can absorb neutrons to become Pu-239)
    % Dark blue circles = U-235 (active fissile material)

    rng(42);  % Fixed seed for reproducibility
    u235_fraction = 0.15;     % 15% of fuel cells contain U-235
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

    if animate
        % Include neutron in legend when animating
        h_neutron_legend = scatter(NaN, NaN, 25, 'k', 'filled');
        legend([h_boundary, h_coolant, h_moderator, h_u235, h_u238, h_neutron_legend], ...
               {sprintf('Fuel Rod Boundary (%d rods)', n_fuel_rods), ...
                'Light Water Coolant', ...
                'Graphite Moderator', ...
                sprintf('U-235 Fissile (%d)', n_u235), ...
                sprintf('U-238 Fertile (%d)', n_fuel_cells - n_u235), ...
                'Neutrons'}, ...
               'Location', 'northeastoutside', ...
               'FontSize', 10);
    else
        legend([h_boundary, h_coolant, h_moderator, h_u235, h_u238], ...
               {sprintf('Fuel Rod Boundary (%d rods)', n_fuel_rods), ...
                'Light Water Coolant', ...
                'Graphite Moderator', ...
                sprintf('U-235 Fissile (%d)', n_u235), ...
                sprintf('U-238 Fertile (%d)', n_fuel_cells - n_u235)}, ...
               'Location', 'northeastoutside', ...
               'FontSize', 10);
    end

    % ========== NEUTRON CHAIN REACTION ANIMATION ==========
    if animate
        fprintf('\nStarting neutron chain reaction animation...\n');

        % Animation parameters
        dt = 0.05;                % Time step per frame
        pause_time = 0.05;        % Pause between frames (slower = more visible)
        neutron_speed = 2.0;      % Grid units per second (slow for visibility)
        collision_radius = 0.35;  % Same as pellet_radius
        max_neutrons = 500;       % Cap to prevent runaway
        n_frames = 800;           % Total animation frames

        % Neutron state arrays
        neutron_x = [];
        neutron_y = [];
        neutron_vx = [];
        neutron_vy = [];

        % Create scatter plot for neutrons (efficient updates)
        h_neutrons = scatter([], [], 30, 'k', 'filled');

        % Start with a single neutron at the left edge, aimed at a random U-235 cell
        target_idx = randi(n_u235);             % Pick a random U-235 as target
        target_x = u235_cx(target_idx);
        target_y = u235_cy(target_idx);

        neutron_x = 0;                          % Start at left edge
        neutron_y = target_y;                   % Align with target U-235 vertically

        % Calculate direction toward target
        dx = target_x - neutron_x;
        dy = target_y - neutron_y;
        dist_to_target = sqrt(dx^2 + dy^2);

        neutron_vx = neutron_speed * (dx / dist_to_target);
        neutron_vy = neutron_speed * (dy / dist_to_target);

        fprintf('Initial neutron entering from left edge, aimed at U-235 at (%.1f, %.1f)\n', target_x, target_y);

        % Main animation loop
        for frame = 1:n_frames
            % Step 1: Move neutrons (smooth sub-step movement)
            neutron_x = neutron_x + neutron_vx * dt;
            neutron_y = neutron_y + neutron_vy * dt;

            % Step 2: Check collisions and handle fission
            new_x = []; new_y = []; new_vx = []; new_vy = [];
            remove = false(size(neutron_x));

            for i = 1:length(neutron_x)
                % Check if out of bounds - remove neutron
                if neutron_x(i) < 0 || neutron_x(i) > n_cols || ...
                   neutron_y(i) < 0 || neutron_y(i) > n_rows
                    remove(i) = true;
                    continue;
                end

                % Check collision with U-235 cells (vectorized distance)
                dist = sqrt((u235_cx - neutron_x(i)).^2 + (u235_cy - neutron_y(i)).^2);
                hit_idx = find(dist < collision_radius, 1);

                if ~isempty(hit_idx) && (length(neutron_x) + length(new_x)) < max_neutrons
                    % Fission: spawn exactly 3 new neutrons
                    [nx, ny, nvx, nvy] = spawn_neutrons(u235_cx(hit_idx), u235_cy(hit_idx), neutron_speed, 3);
                    new_x = [new_x; nx];
                    new_y = [new_y; ny];
                    new_vx = [new_vx; nvx];
                    new_vy = [new_vy; nvy];
                    remove(i) = true;  % Remove triggering neutron
                end
            end

            % Remove collided/escaped neutrons and add new ones
            neutron_x = [neutron_x(~remove); new_x];
            neutron_y = [neutron_y(~remove); new_y];
            neutron_vx = [neutron_vx(~remove); new_vx];
            neutron_vy = [neutron_vy(~remove); new_vy];

            % Step 3: Update graphics smoothly
            set(h_neutrons, 'XData', neutron_x, 'YData', neutron_y);
            drawnow;
            pause(pause_time);

            % End animation if no neutrons left - restart from left edge aimed at U-235
            if isempty(neutron_x)
                fprintf('All neutrons escaped. Restarting with new neutron...\n');
                target_idx = randi(n_u235);
                target_x = u235_cx(target_idx);
                target_y = u235_cy(target_idx);

                neutron_x = 0;
                neutron_y = target_y;

                dx = target_x - neutron_x;
                dy = target_y - neutron_y;
                dist_to_target = sqrt(dx^2 + dy^2);

                neutron_vx = neutron_speed * (dx / dist_to_target);
                neutron_vy = neutron_speed * (dy / dist_to_target);
            end

            % Progress indicator every 200 frames
            if mod(frame, 200) == 0
                fprintf('Frame %d/%d - Active neutrons: %d\n', frame, n_frames, length(neutron_x));
            end
        end

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
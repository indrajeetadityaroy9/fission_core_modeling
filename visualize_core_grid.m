function fig = visualize_core_grid(varargin)
%VISUALIZE_CORE_GRID Creates a grid cross-section of the RBMK reactor core
%
%   fig = visualize_core_grid()
%   fig = visualize_core_grid('preset', 'high_power')
%
%   Displays a cross-section view of the reactor core as a grid.
%   Interleaving vertical rectangles of fuel rods and graphite moderators.
%   Light blue squares represent coolant in fuel rod channels.
%
%   Optional Parameters:
%       'preset'   - Configuration preset: 'high_power' (default: 'high_power')
%       'rods_out' - If true, all control rods are lifted (no absorption),
%                    disabling automatic control (default: false)
%
%   Output:
%       fig - Figure handle

    % Parse optional parameters
    parser = inputParser;
    addParameter(parser, 'preset', 'high_power', @ischar);
    addParameter(parser, 'rods_out', false, @islogical);  % Lift all control rods (no absorption)
    parse(parser, varargin{:});
    preset   = parser.Results.preset;
    rods_out = parser.Results.rods_out;

    % ========== LOAD CONFIGURATION PRESET ==========
    cfg = get_preset_config(preset);

    % ========== DERIVED VALUES ==========
    n_rows             = cfg.grid.n_rows;
    n_fuel_strips      = cfg.grid.n_fuel_strips;
    n_moderator_strips = cfg.grid.n_moderator_strips;
    n_cols = n_fuel_strips * cfg.grid.fuel_width + ...
             n_moderator_strips * cfg.grid.moderator_width;
    moderator_width = cfg.grid.moderator_width;
    fuel_width      = cfg.grid.fuel_width;

    % Create a map: true = fuel rod (coolant), false = moderator
    fuel_map = false(n_rows, n_cols);

    % Build interleaving pattern: M(3) | F(5) | M(3) | F(5) | ...
    col = 1;
    fuel_strip_positions = [];       % Track fuel strip start positions for boundaries
    moderator_strip_positions = [];  % Track moderator strip start positions for rods
    for i = 1:(n_fuel_strips + n_moderator_strips)
        if mod(i, 2) == 1
            % Odd = moderator strip
            moderator_strip_positions = [moderator_strip_positions; col];
            col = col + moderator_width;
        else
            % Even = fuel strip
            fuel_strip_positions = [fuel_strip_positions; col];
            fuel_map(:, col:col+fuel_width-1) = true;
            col = col + fuel_width;
        end
    end

    n_fuel_rods = cfg.grid.n_fuel_strips;

    % Control rod positions (center column of each moderator strip)
    n_control_rods   = numel(moderator_strip_positions);
    rod_center_offset = floor((moderator_width - 1) / 2);
    control_rod_cols = moderator_strip_positions - 0.5 + rod_center_offset;
    control_rod_width = cfg.fuel.control_rod_width;

    % Rod position state (0 = fully raised, n_rows = fully inserted)
    % Rods enter from top (y=0) and extend downward
    if rods_out
        rod_position = zeros(n_control_rods, 1);    % All rods fully raised
        rod_target   = zeros(n_control_rods, 1);    % Target: stay raised
    else
        rod_position = n_rows * ones(n_control_rods, 1);  % All rods fully inserted
        rod_target   = n_rows * ones(n_control_rods, 1);  % Target: stay inserted
    end

    % Coolant state arrays (for fuel channel cells where fuel_map is true)
    coolant_temp   = zeros(n_rows, n_cols);      % Temperature (0 = cold, 1 = boiling)
    coolant_density = ones(n_rows, n_cols);      % Density (1 = liquid, 0 = void/steam)
    coolant_void   = false(n_rows, n_cols);      % True when evaporated

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
                cell_color = cfg.colors.coolant;
            else
                cell_color = cfg.colors.moderator;
            end

            cell_patches(row, col) = patch(x, y, cell_color, ...
                  'EdgeColor', cfg.colors.edge, ...
                  'LineWidth', 0.3);
        end
    end

    % Draw control rods in moderator strips (boron absorbers)
    % Height based on rod_position (0 = invisible, n_rows = full insertion)
    h_control_rods = gobjects(n_control_rods, 1);  % Graphics handles for dynamic updates
    for i = 1:n_control_rods
        x_left = control_rod_cols(i) - control_rod_width/2;
        h_control_rods(i) = rectangle('Position', [x_left, 0, control_rod_width, rod_position(i)], ...
                                       'FaceColor', cfg.colors.control_rod, ...
                                       'EdgeColor', 'k', ...
                                       'LineWidth', 1);
    end

    % Add uranium circles to all fuel rod cells
    % Grey circles = U-238 (fertile material, can absorb neutrons to become Pu-239)
    % Yellow circles = U-235 (active fissile material)

    rng(42);  % Fixed seed for reproducibility
    u235_fraction  = cfg.fuel.u235_fraction;
    pellet_radius  = cfg.fuel.pellet_radius;

    % Get all fuel cell positions
    [fuel_rows, fuel_cols] = find(fuel_map);
    n_fuel_cells = length(fuel_rows);
    n_u235       = round(u235_fraction * n_fuel_cells);

    % Randomly select cells for U-235
    u235_indices = randperm(n_fuel_cells, n_u235);
    u235_set     = false(n_fuel_cells, 1);
    u235_set(u235_indices) = true;

    % Draw circles for all fuel cells
    % Store U-235 patch handles and ring handles for xenon poisoning visualization
    u235_patches     = gobjects(n_u235, 1);
    u235_rings       = gobjects(n_u235, 1);   % Ring indicators for xenon poisoning
    u235_xe_overlays = gobjects(n_u235, 1);   % Xenon overlay patches
    u235_patch_idx = 1;

    theta       = linspace(0, 2*pi, 30);
    ring_radius = pellet_radius + 0.08;  % Ring slightly larger than pellet

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
            u235_patches(u235_patch_idx) = patch(rx, ry, cfg.colors.uranium, ...
                  'EdgeColor', cfg.colors.uranium * 0.7, ...
                  'LineWidth', 0.5);

            % Xenon overlay: same footprint, starts transparent
            u235_xe_overlays(u235_patch_idx) = patch(rx, ry, cfg.colors.xenon, ...
                  'EdgeColor', 'none', ...
                  'FaceAlpha', 0);

            % Draw xenon poisoning ring indicator (initially hidden)
            ring_x = cx + ring_radius * cos(theta);
            ring_y = cy + ring_radius * sin(theta);
            u235_rings(u235_patch_idx) = patch(ring_x, ring_y, [1 1 1], ...
                  'FaceColor', 'none', ...
                  'EdgeColor', cfg.colors.xenon, ...
                  'LineWidth', 2.5, ...
                  'Visible', 'off');

            u235_patch_idx = u235_patch_idx + 1;
        else
            % U-238 (fertile material) - grey
            patch(rx, ry, cfg.colors.u238, ...
                  'EdgeColor', cfg.colors.u238 * 0.7, ...
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
    % KD-tree for fast nearest-neighbor queries against U-235 sites
    u235_tree = createns([u235_cx, u235_cy], 'NSMethod', 'kdtree');

    % Upper/lower region split (axial)
    mid_row    = n_rows / 2;
    upper_rows = 1:mid_row;
    lower_rows = mid_row+1:n_rows;

    % Optional: classify U-235 sites by region
    idx_upper_u235 = find(u235_cy <= mid_row);
    idx_lower_u235 = find(u235_cy >  mid_row);

    % Xenon-135 poisoning state (per U-235 site)
    xenon_level = zeros(n_u235, 1);  % Poison level (0 = clean, 1 = fully poisoned)
    u235_active = true(n_u235, 1);   % false = poisoned (no longer fissile)

    % Configure axes
    xlim([0, n_cols]);
    ylim([0, n_rows]);

    % Flip y-axis so row 1 is at top
    set(gca, 'YDir', 'reverse');

    % Labels
    xlabel('Column', 'FontSize', 12);
    ylabel('Row',    'FontSize', 12);
    title(sprintf('RBMK-1000 Reactor Core Cross-Section (%d × %d Grid)', n_rows, n_cols), ...
          'FontSize', 14, 'FontWeight', 'bold');

    % Add tick marks
    set(gca, 'XTick', 0:5:n_cols, 'YTick', 0:6:n_rows);
    set(gca, 'TickDir', 'out');
    set(gca, 'FontSize', 10);

    % Visual separator between upper and lower regions
    sep_y = mid_row;
    plot([0, n_cols], [sep_y, sep_y], '--', ...
         'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    text(0.3, mid_row/2, ...
         'Upper region', 'Color','k', 'FontSize',11, 'FontWeight','bold');
    text(0.3, mid_row + (n_rows - mid_row)/2, ...
         'Lower region', 'Color','k', 'FontSize',11, 'FontWeight','bold');

    % ========== NEUTRON CHAIN REACTION ANIMATION ==========
    fprintf('\nStarting neutron chain reaction animation...\n');

    % Extract animation parameters from cfg
    dt             = cfg.anim.dt;
    pause_time     = cfg.anim.pause_time;
    neutron_speed  = cfg.anim.neutron_speed;
    collision_radius = cfg.anim.collision_radius;
    max_neutrons   = cfg.anim.max_neutrons;
    n_frames       = cfg.anim.n_frames;
    rod_speed      = cfg.anim.rod_speed;
    substeps       = cfg.anim.substeps;
    sub_dt         = dt / substeps;   % Sub-stepping to reduce tunneling

    % Extract physics parameters from cfg
    neutron_threshold       = cfg.physics.neutron_threshold;
    temp_baseline           = cfg.physics.temp_baseline;
    temp_evaporate          = cfg.physics.temp_evaporate;
    temp_condense           = cfg.physics.temp_condense;
    heat_per_neutron        = cfg.physics.heat_per_neutron;
    cooling_rate            = cfg.physics.cooling_rate;
    max_coolant_absorption  = cfg.physics.max_coolant_absorption;
    xenon_per_fission       = cfg.physics.xenon_per_fission;
    xenon_threshold         = cfg.physics.xenon_threshold;
    xenon_decay_rate        = cfg.physics.xenon_decay_rate;
    spontaneous_fission_prob = cfg.physics.spontaneous_fission_prob;

    % Neutron balance tracking
    total_produced        = 1;   % Start with 1 initial neutron
    total_absorbed_rods   = 0;
    total_absorbed_coolant = 0;
    total_absorbed_xenon  = 0;
    total_fissions        = 0;
    total_escaped         = 0;
    total_spontaneous     = 0;

    % Neutron state arrays
    neutron_x  = [];
    neutron_y  = [];
    neutron_vx = [];
    neutron_vy = [];
    neutron_in_rod = [];  % Track current control-rod index a neutron is inside (0 = none)

    % Create scatter plot for neutrons
    h_neutrons = scatter([], [], 30, 'k', 'filled');

    % Start with a single neutron aimed at a random U-235 cell
    target_idx = randi(n_u235);             % Pick a random U-235 as target
    target_x   = u235_cx(target_idx);
    target_y   = u235_cy(target_idx);

    % Start just inside the first fuel rod (past the left moderator strip)
    neutron_x = [moderator_width + 0.5];    % Start inside first fuel rod
    neutron_y = [target_y];                 % Align with target U-235 vertically
    neutron_in_rod = zeros(size(neutron_x));

    % Calculate direction toward target
    dx = target_x - neutron_x;
    dy = target_y - neutron_y;
    dist_to_target = hypot(dx, dy);

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
        % Step 0: Automatic control logic - set target positions for control rods
        % (disabled when rods_out mode is active - all rods stay lifted)
        if ~rods_out
            active_neutron_count = length(neutron_x);

            if active_neutron_count < neutron_threshold
                % Set target to RAISED for alternating rods (2 and 4) to allow more fission
                for j = 2:2:n_control_rods
                    if rod_target(j) ~= 0
                        rod_target(j) = 0;  % Target: fully raised
                        fprintf('Rod %d RAISING (neutrons: %d < %d)\n', j, active_neutron_count, neutron_threshold);
                    end
                end
            elseif active_neutron_count > neutron_threshold
                % Set target to INSERTED for all rods to suppress chain reaction
                for j = 1:n_control_rods
                    if rod_target(j) ~= n_rows
                        rod_target(j) = n_rows;  % Target: fully inserted
                        fprintf('Rod %d INSERTING (neutrons: %d > %d)\n', j, active_neutron_count, neutron_threshold);
                    end
                end
            end
        end

        % Step 0.25: Animate control rod movement toward target positions
        for j = 1:n_control_rods
            if rod_position(j) < rod_target(j)
                % Inserting (moving down into core)
                rod_position(j) = min(rod_target(j), rod_position(j) + rod_speed);
            elseif rod_position(j) > rod_target(j)
                % Raising (withdrawing from core)
                rod_position(j) = max(rod_target(j), rod_position(j) - rod_speed);
            end

            % Update rectangle height (grows from top y=0 downward)
            x_left = control_rod_cols(j) - control_rod_width/2;
            set(h_control_rods(j), 'Position', [x_left, 0, control_rod_width, max(0.01, rod_position(j))]);
            end

            % Step 0.5: Background cooling for all coolant cells (once per frame)
            for row = 1:n_rows
                for col = 1:n_cols
                if fuel_map(row, col)
                    % Gradual cooling toward baseline
                    coolant_temp(row, col) = max(temp_baseline, ...
                        coolant_temp(row, col) - cooling_rate);

                    % Check for condensation (void → liquid)
                    if coolant_void(row, col) && coolant_temp(row, col) < temp_condense
                        coolant_void(row, col)   = false;
                        coolant_density(row, col) = 1.0;
                    end

                    % Update density based on temperature (gradual reduction as heating)
                        if ~coolant_void(row, col)
                            coolant_density(row, col) = max(0, 1 - coolant_temp(row, col));
                        end
                    end
                end
            end

            % Step 0.6: Void convection (bubbles drift upward)
            void_rise_prob = 0.7;
            for col = 1:n_cols
                for row = n_rows:-1:2  % bottom to top
                    if ~fuel_map(row, col)
                        continue;
                    end
                    if coolant_void(row, col) && ~coolant_void(row-1, col)
                        if rand() < void_rise_prob
                            % Swap states with cell above
                            tmp_void       = coolant_void(row-1, col);
                            tmp_temp       = coolant_temp(row-1, col);
                            tmp_density    = coolant_density(row-1, col);

                            coolant_void(row-1, col)    = coolant_void(row, col);
                            coolant_temp(row-1, col)    = coolant_temp(row, col);
                            coolant_density(row-1, col) = coolant_density(row, col);

                            coolant_void(row, col)    = tmp_void;
                            coolant_temp(row, col)    = tmp_temp;
                            coolant_density(row, col) = tmp_density;
                        end
                    end
                end
            end

            % Remove voids that reach the very top (leave the core)
            for col = 1:n_cols
                if fuel_map(1, col) && coolant_void(1, col)
                    coolant_void(1, col)    = false;
                    coolant_density(1, col) = 1.0;
                    coolant_temp(1, col)    = temp_baseline;
                end
            end

            % ===== SUB-STEPPED NEUTRON MOTION & COLLISIONS =====
            for s = 1:substeps
            % Step 1: Move neutrons (smaller sub-step)
            neutron_x = neutron_x + neutron_vx * sub_dt;
            neutron_y = neutron_y + neutron_vy * sub_dt;

                % Step 2: Check collisions and handle physics for this sub-step
                new_x = []; new_y = []; new_vx = []; new_vy = [];
                remove = false(size(neutron_x));
                if isempty(neutron_x)
                    u235_neighbors = cell(0, 1);
                else
                    u235_neighbors = rangesearch(u235_tree, [neutron_x(:) neutron_y(:)], collision_radius);
                end

                for i = 1:length(neutron_x)
                    % Check if out of bounds - remove neutron (escaped)
                    if neutron_x(i) < 0 || neutron_x(i) > n_cols || ...
                       neutron_y(i) < 0 || neutron_y(i) > n_rows
                    remove(i) = true;
                    total_escaped = total_escaped + 1;
                    continue;
                end

                % Check collision with control rods (only in inserted portion)
                % Control rod interaction: apply absorption once per rod crossing
                current_rod = 0;
                for j = 1:n_control_rods
                    if rod_position(j) > 0
                        rod_x        = control_rod_cols(j);
                        rod_half_width = control_rod_width / 2;
                        if neutron_x(i) >= rod_x - rod_half_width && ...
                           neutron_x(i) <= rod_x + rod_half_width && ...
                           neutron_y(i) <= rod_position(j)
                            current_rod = j;
                            break;
                        end
                    end
                end

                if current_rod > 0 && neutron_in_rod(i) ~= current_rod
                    depth_fraction  = min(1, rod_position(current_rod) / n_rows);
                    absorb_prob     = cfg.physics.rod_absorption_max * depth_fraction;
                    if rand() < absorb_prob
                        remove(i) = true;
                        total_absorbed_rods = total_absorbed_rods + 1;
                        continue;
                    else
                        neutron_in_rod(i) = current_rod;  % Mark as traversing this rod
                    end
                elseif current_rod == 0
                    neutron_in_rod(i) = 0;  % Reset when outside any rod
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
                coolant_void(cell_row, cell_col)   = true;
                coolant_density(cell_row, cell_col) = 0;

                % High-power Doppler brake: voids are swept out immediately
                if cfg.physics.clear_void_immediately
                    coolant_void(cell_row, cell_col)   = false;
                    coolant_density(cell_row, cell_col) = 1.0;
                    coolant_temp(cell_row, cell_col)    = temp_baseline;
                end
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
                    hit_candidates = u235_neighbors{i};
                    if ~isempty(hit_candidates)
                        hit_idx = hit_candidates(1);  % Take first nearby site
                        if ~u235_active(hit_idx)
                            % Poisoned site (xenon-135) - absorb neutron and burn off xenon
                            remove(i) = true;
                            total_absorbed_xenon = total_absorbed_xenon + 1;

                        % Neutron absorption burns off xenon (Xe-135 + n -> Xe-136)
                        xenon_level(hit_idx) = xenon_level(hit_idx) - xenon_per_fission * 2;

                        % Check if site recovers from poisoning
                        if xenon_level(hit_idx) < xenon_threshold * 0.5
                            u235_active(hit_idx) = true;
                            xenon_level(hit_idx) = max(0, xenon_level(hit_idx));
                            set(u235_rings(hit_idx), 'Visible', 'off');  % Hide ring
                        end
                    else
                        % Active site - fission if under neutron limit
                        current_count = length(neutron_x) - sum(remove) + length(new_x);
                        if current_count < max_neutrons
                            % Fission: spawn 3 neutrons
                            [nx, ny, nvx, nvy] = spawn_neutrons(u235_cx(hit_idx), u235_cy(hit_idx), neutron_speed, 3);
                            new_x  = [new_x;  nx];
                            new_y  = [new_y;  ny];
                            new_vx = [new_vx; nvx];
                            new_vy = [new_vy; nvy];
                            remove(i) = true;  % Remove triggering neutron
                            total_fissions  = total_fissions + 1;
                            total_produced  = total_produced + 3;

                            % Accumulate xenon from fission
                            xenon_level(hit_idx) = xenon_level(hit_idx) + xenon_per_fission;

                            % Check if site becomes poisoned
                            if xenon_level(hit_idx) >= xenon_threshold
                                u235_active(hit_idx) = false;
                            end
                        end
                    end
                end
            end

                % Remove absorbed/escaped neutrons and add new ones
                keep_idx        = find(~remove);
                neutron_in_rod  = [neutron_in_rod(keep_idx); zeros(numel(new_x), 1)];
                neutron_x       = [neutron_x(keep_idx);  new_x(:)];
                neutron_y       = [neutron_y(keep_idx);  new_y(:)];
                neutron_vx      = [neutron_vx(keep_idx); new_vx(:)];
                neutron_vy      = [neutron_vy(keep_idx); new_vy(:)];
            end
        % ===== END SUB-STEPS =====

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

        % Step 3.5: Update U-235 colors for xenon poisoning
        for i = 1:n_u235
            % Optional: slow xenon decay (burnout) for active sites only
            if xenon_level(i) > 0 && u235_active(i)
                xenon_level(i) = max(0, xenon_level(i) - xenon_decay_rate);
            end

            % Xenon overlay opacity based on poisoning (0 = clean, 1 = fully poisoned)
            poison_frac = min(1.0, xenon_level(i) / xenon_threshold);
            alpha_max   = 0.85;
            xe_alpha    = poison_frac * alpha_max;
            set(u235_xe_overlays(i), 'FaceAlpha', xe_alpha);

            % Rings remain as a binary indicator for fully poisoned sites
            if ~u235_active(i)
                set(u235_rings(i), 'Visible', 'on');
            else
                set(u235_rings(i), 'Visible', 'off');
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

        % Spontaneous fission - background neutron source from fissile material
        if rand() < spontaneous_fission_prob && length(neutron_x) < max_neutrons
            % Find active (non-poisoned) U-235 sites for spontaneous fission
            active_sites = find(u235_active);
            if ~isempty(active_sites)
                site_idx = active_sites(randi(length(active_sites)));
                [sx, sy, svx, svy] = spawn_neutrons(u235_cx(site_idx), u235_cy(site_idx), neutron_speed, 1);
                neutron_x  = [neutron_x;  sx];
                neutron_y  = [neutron_y;  sy];
                neutron_vx = [neutron_vx; svx];
                neutron_vy = [neutron_vy; svy];
                total_spontaneous = total_spontaneous + 1;
                total_produced    = total_produced + 1;
            end
        end

        % Progress indicator every 100 frames with balance info
        if mod(frame, 100) == 0
            n_voids    = sum(coolant_void(:));
            n_poisoned = sum(~u235_active);

            voids_upper = sum(sum(coolant_void(upper_rows, :) & fuel_map(upper_rows, :)));
            voids_lower = sum(sum(coolant_void(lower_rows, :) & fuel_map(lower_rows, :)));

            fprintf(['Frame %d/%d | Active: %d | Fissions: %d | ', ...
                     'Voids=%d (U=%d, L=%d) | Poisoned: %d/%d\n'], ...
                    frame, n_frames, length(neutron_x), ...
                    total_fissions, ...
                    n_voids, voids_upper, voids_lower, ...
                    n_poisoned, n_u235);
        end
    end

    % Final balance summary
    n_voids    = sum(coolant_void(:));
    n_poisoned = sum(~u235_active);
    fprintf('\n===== NEUTRON BALANCE SUMMARY =====\n');
    fprintf('Total produced:      %d\n', total_produced);
    fprintf('  - From fission:    %d\n', total_fissions * 3);
    fprintf('  - Spontaneous:     %d\n', total_spontaneous);
    fprintf('Absorbed (coolant):  %d\n', total_absorbed_coolant);
    fprintf('Absorbed (rods):     %d\n', total_absorbed_rods);
    fprintf('Absorbed (xenon):    %d\n', total_absorbed_xenon);
    fprintf('Total escaped:       %d (boundary)\n', total_escaped);
    fprintf('Total fissions:      %d\n', total_fissions);
    fprintf('Active remaining:    %d\n', length(neutron_x));
    fprintf('Void cells:          %d / %d (%.1f%%)\n', n_voids, n_fuel_cells, 100*n_voids/n_fuel_cells);
    fprintf('Poisoned U-235:      %d / %d (%.1f%%)\n', n_poisoned, n_u235, 100*n_poisoned/n_u235);
    fprintf('Neutron limit:       %d (max allowed)\n', max_neutrons);
    fprintf('===================================\n');
    fprintf('Animation complete.\n');

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

function cfg = get_preset_config(preset_name)
%GET_PRESET_CONFIG Returns configuration struct for the specified preset
%   preset_name - Name of preset: 'high_power'

    switch preset_name
        case 'high_power'
            cfg = high_power_config();
        otherwise
            warning('Unknown preset "%s", using high_power', preset_name);
            cfg = high_power_config();
    end
end

function cfg = high_power_config()
%HIGH_POWER_CONFIG Normal high-power operation (~3000 MW thermal)
%   Dense neutron cloud, strong cooling flow, stable operation
%   Simulates reactor functioning as designed with LAC system active

    cfg = struct();

    % --- Grid dimensions (Standard) ---
    cfg.grid.n_rows             = 24;
    cfg.grid.moderator_width    = 3;
    cfg.grid.fuel_width         = 5;
    cfg.grid.n_fuel_strips      = 4;
    cfg.grid.n_moderator_strips = 5;

    % --- Colors (Standard) ---
    cfg.colors.coolant      = [0.15 0.25 0.7];
    cfg.colors.moderator    = [0.85 0.85 0.85];
    cfg.colors.edge         = [0.5 0.5 0.5];
    cfg.colors.fuel_boundary = [0.8 0.2 0.1];
    cfg.colors.uranium      = [1.0 0.88 0.2];  % U-235 pellets (yellow)
    cfg.colors.u238         = [0.5 0.5 0.5];
    cfg.colors.control_rod  = [0.2 0.2 0.2];
    cfg.colors.xenon        = [1.0 0.5 0.0];

    % --- Fuel configuration ---
    cfg.fuel.u235_fraction     = 0.35;
    cfg.fuel.pellet_radius     = 0.35;
    cfg.fuel.control_rod_width = 0.8;

    % --- Animation parameters (Optimized for High Activity) ---
    cfg.anim.dt             = 0.05;
    cfg.anim.pause_time     = 0.03;    % Faster animation speed
    cfg.anim.neutron_speed  = 3.0;     % Faster neutrons for high flux
    cfg.anim.collision_radius = 0.35;
    cfg.anim.max_neutrons   = 300;     % HIGH LIMIT: Dense neutron cloud
    cfg.anim.n_frames       = 1000;
    cfg.anim.rod_speed      = 0.8;     % Faster rod movement for high-activity sim
    cfg.anim.substeps       = 4;       % Substeps for motion to reduce tunneling

    % --- Physics parameters (STABLE HIGH-POWER REGIME) ---
    % High threshold simulates operating at ~3000 MW thermal
    cfg.physics.neutron_threshold = 180;

    % Doppler + fast-flow feedback: voids are cleared immediately at high power
    cfg.physics.clear_void_immediately = true;

    % Control rod absorption (linear with insertion depth)
    cfg.physics.rod_absorption_max = 0.9;

    % Strong cooling capacity - heat removed quickly, prevents void formation
    cfg.physics.cooling_rate = 0.009;

    % Standard heat generation per neutron
    cfg.physics.heat_per_neutron = 0.03;

    % Thermal hydraulics - high pressure (7 MPa) keeps water liquid
    cfg.physics.temp_baseline = 0.0;
    cfg.physics.temp_evaporate = 0.95;
    cfg.physics.temp_condense  = 0.25;

    % Water absorption provides negative reactivity feedback
    cfg.physics.max_coolant_absorption = 0.008;  % Low per-substep rate

    % Xenon dynamics - high power burns off xenon (equilibrium)
    cfg.physics.xenon_per_fission = 0.15;
    cfg.physics.xenon_threshold   = 1.0;
    cfg.physics.xenon_decay_rate  = 0.005;  % Fast decay (burnoff)

    % Spontaneous fission - background neutron source
    cfg.physics.spontaneous_fission_prob = 0.03;  % Slightly higher at high power

    % Void handling: allow visualization of upward drift (do not clear immediately)
    cfg.physics.clear_void_immediately = false;
end
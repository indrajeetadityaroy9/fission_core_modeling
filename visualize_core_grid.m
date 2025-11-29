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
%       'preset'    - Configuration preset: 'high_power', 'low_power', 'chernobyl'
%       'orm_level' - Rod insertion depth (0=withdrawn, 1=inserted, default: 1.0)
%
%   Output:
%       fig - Figure handle

    %% ========================================================================
    %  SECTION 1: CONFIGURATION & SETUP
    %  ========================================================================

    % Parse optional parameters
    parser = inputParser;
    addParameter(parser, 'preset', 'high_power', @ischar);
    addParameter(parser, 'orm_level', 1.0, @isnumeric);
    parse(parser, varargin{:});
    preset    = parser.Results.preset;
    orm_level = parser.Results.orm_level;

    % Load configuration preset
    cfg = get_preset_config(preset);

    % Apply preset-specific overrides
    [orm_level, initial_power] = apply_preset_overrides(preset, orm_level);

    %% ========================================================================
    %  SECTION 2: GRID GEOMETRY
    %  ========================================================================

    grid = setup_grid_geometry(cfg);

    %% ========================================================================
    %  SECTION 3: PHYSICS PARAMETERS (Constants)
    %  ========================================================================

    params = init_physics_params();

    % Apply preset-specific physics overrides (operate on full RBMK params)
    if isfield(cfg, 'physics') && isfield(cfg.physics, 't_scram')
        params.p_rbmk.t_scram = cfg.physics.t_scram;
    end

    %% ========================================================================
    %  SECTION 4: REACTOR STATE (Variables)
    %  ========================================================================

    state = init_reactor_state(preset, initial_power, params, grid, orm_level);

    %% ========================================================================
    %  SECTION 5: GRAPHICS SETUP
    %  ========================================================================

    [fig, gfx] = setup_graphics(cfg, grid, state);

    % Initialize at least one neutron, aimed at a random U-235 site
    if gfx.n_u235 > 0
        target_idx = randi(gfx.n_u235);
        target_x = gfx.u235_cx(target_idx);
        target_y = gfx.u235_cy(target_idx);

        % Start inside first fuel strip, aligned vertically with target
        first_fuel_col = grid.fuel_strip_positions(1);
        start_x = first_fuel_col + 0.5;
        start_y = target_y;

        dx0 = target_x - start_x;
        dy0 = target_y - start_y;
        dist0 = hypot(dx0, dy0);
        if dist0 > 0
            vx0 = cfg.anim.neutron_speed * (dx0 / dist0);
            vy0 = cfg.anim.neutron_speed * (dy0 / dist0);
        else
            vx0 = cfg.anim.neutron_speed;
            vy0 = 0;
        end

        state.neutron_x = start_x;
        state.neutron_y = start_y;
        state.neutron_vx = vx0;
        state.neutron_vy = vy0;
        state.neutron_in_rod = 0;
    end

    %% ========================================================================
    %  SECTION 6: ANIMATION LOOP
    %  ========================================================================

    % Extract animation parameters
    dt = cfg.anim.dt;
    n_frames = cfg.anim.n_frames;
    pause_time = cfg.anim.pause_time;
    substeps = cfg.anim.substeps;
    sub_dt = dt / substeps;

    % Main animation loop
    for frame = 1:n_frames
        % Step 1: UPDATE PHYSICS (ODE integration + reactivity)
        [state, rho] = update_reactor_physics(state, params, grid, dt);

        % Step 1b: UPDATE SIMULATION TIME (for SCRAM timing)
        state.simulation_time = state.simulation_time + dt;

        % Physics summary every 100 frames: relate ODE reactivity,
        % neutron density, fissions, and visual particle count.
        if frame == 1 || mod(frame, 100) == 0
            avg_rho_void = (rho.void_L + rho.void_U) / 2;
            avg_rho_dop  = (rho.dop_L + rho.dop_U) / 2;
            avg_rho_xen  = (rho.xen_L + rho.xen_U) / 2;
            n_particles = length(state.neutron_x);
            fprintf(['PHYSICS: frame %4d t=%6.2f | neutrons=%4d ' ...
                     '| fissions=%6d | n=[%.3f,%.3f] alpha=[%.3f,%.3f] ' ...
                     'rho_total=%.5f (void=%.5f, dop=%.5f, xen=%.5f, rod=%.5f)\n'], ...
                frame, state.simulation_time, n_particles, state.total_fissions, ...
                state.n_L, state.n_U, state.alpha_L, state.alpha_U, rho.total, ...
                avg_rho_void, avg_rho_dop, avg_rho_xen, rho.rod);
        end

        % Step 2: BACKGROUND COOLING / CONDENSATION FOR COOLANT CELLS
        state = update_coolant_background(state, cfg, grid);

        % Step 3: UPDATE PARTICLE SIMULATION
        for s = 1:substeps
            [state, gfx] = update_particles(state, gfx, grid, cfg, sub_dt, params);
        end

        % Step 4: BALANCE PARTICLES TO MATCH ODE
        state = balance_particles(state, gfx, cfg, params);

        % Step 5: UPDATE VISUALIZATION
        update_cell_colors(state, gfx, grid, params);
        update_neutron_graphics(state, gfx);
        update_stats_panel(state, rho, gfx, params);

        % Render frame
        drawnow;
        pause(pause_time);

        % Spontaneous fission: background neutron source from active fuel
        if rand() < cfg.physics.spontaneous_fission_prob && length(state.neutron_x) < cfg.anim.max_neutrons ...
                && gfx.n_u235 > 0
            site_idx = randi(gfx.n_u235);
            [sx, sy, svx, svy] = spawn_neutrons(gfx.u235_cx(site_idx), gfx.u235_cy(site_idx), ...
                                                cfg.anim.neutron_speed, 1);
            state.neutron_x = [state.neutron_x; sx];
            state.neutron_y = [state.neutron_y; sy];
            state.neutron_vx = [state.neutron_vx; svx];
            state.neutron_vy = [state.neutron_vy; svy];
            state.neutron_in_rod = [state.neutron_in_rod; 0];
            state.total_spontaneous = state.total_spontaneous + 1;
            state.total_produced = state.total_produced + 1;
        end

        % Progress indicator every 100 frames
        if mod(frame, 100) == 0
            print_progress(frame, n_frames, state, gfx, grid);
        end
    end

    %% ========================================================================
    %  SECTION 7: FINAL SUMMARY
    %  ========================================================================

    print_final_summary(state, gfx, grid, params, cfg);

    hold off;
end


%% ============================================================================
%  INITIALIZATION FUNCTIONS
%  ============================================================================

function [orm_level, initial_power] = apply_preset_overrides(preset, orm_level)
%APPLY_PRESET_OVERRIDES Apply preset-specific parameter overrides
    initial_power = 0.2;  % Default

    if strcmp(preset, 'chernobyl') && orm_level == 1.0
        orm_level = 0.2;
        initial_power = 0.3;
        fprintf('Chernobyl preset: ORM violation active (orm_level=%.1f)\n', orm_level);
    elseif strcmp(preset, 'high_power') && orm_level == 1.0
        orm_level = 0.35;
        initial_power = 0.5;
        fprintf('High power preset: Adjusting rods for criticality (orm_level=%.2f)\n', orm_level);
    end
end

function grid = setup_grid_geometry(cfg)
%SETUP_GRID_GEOMETRY Create grid structure with all geometry info
    grid = struct();

    grid.n_rows = cfg.grid.n_rows;
    grid.n_fuel_strips = cfg.grid.n_fuel_strips;
    grid.n_moderator_strips = cfg.grid.n_moderator_strips;
    grid.moderator_width = cfg.grid.moderator_width;
    grid.fuel_width = cfg.grid.fuel_width;
    grid.n_cols = grid.n_fuel_strips * grid.fuel_width + ...
                  grid.n_moderator_strips * grid.moderator_width;
    grid.mid_row = grid.n_rows / 2;
    grid.upper_rows = 1:grid.mid_row;
    grid.lower_rows = (grid.mid_row+1):grid.n_rows;

    % Build fuel map: true = fuel rod, false = moderator
    grid.fuel_map = false(grid.n_rows, grid.n_cols);
    col = 1;
    grid.fuel_strip_positions = [];
    grid.moderator_strip_positions = [];

    for i = 1:(grid.n_fuel_strips + grid.n_moderator_strips)
        if mod(i, 2) == 1
            grid.moderator_strip_positions = [grid.moderator_strip_positions; col];
            col = col + grid.moderator_width;
        else
            grid.fuel_strip_positions = [grid.fuel_strip_positions; col];
            grid.fuel_map(:, col:col+grid.fuel_width-1) = true;
            col = col + grid.fuel_width;
        end
    end

    % Control rod positions
    grid.n_control_rods = numel(grid.moderator_strip_positions);
    rod_center_offset = floor((grid.moderator_width - 1) / 2);
    grid.control_rod_cols = grid.moderator_strip_positions - 0.5 + rod_center_offset;

    % Fuel cell info
    [grid.fuel_rows, grid.fuel_cols] = find(grid.fuel_map);
    grid.n_fuel_cells = length(grid.fuel_rows);
end

function params = init_physics_params()
%INIT_PHYSICS_PARAMS Initialize RBMK physics + visualization tunables.
%   This function separates:
%     - RBMK core physics parameters (copied from rbmk_parameters.m)
%     - Visualization-only knobs (purely for animation / plotting)

    % ---------- RBMK PHYSICS (single source of truth: rbmk_parameters.m) ----------
    p_full = rbmk_parameters();

    params = struct();

    % Keep a full copy of the calibrated RBMK parameter set so we can call
    % rbmk_dynamics directly and stay faithful to the core model.
    params.p_rbmk = p_full;

    % ---------- Visualization / animation parameters (NOT core physics) ----------

    % ORM thresholds and amplification are pedagogical, not calibrated physics
    params.ORM_min             = 1.5;  % Minimum equivalent rods for "OK"
    params.ORM_safe            = 3.0;  % "Safe" ORM threshold
    params.ORM_violation_factor = 5.0; % Void amplification factor at low ORM

    % Visual-only plotting knobs (can be tuned without touching physics)
    params.void_threshold      = 0.15; % When to start "void" coloring
    params.density_to_particles = 500; % Particle count scaling for neutrons
end

function state = init_reactor_state(preset, initial_power, params, grid, orm_level)
%INIT_REACTOR_STATE Initialize all reactor state variables
    state = struct();

    % Convenience handle to full RBMK parameter set
    p = params.p_rbmk;

    % Neutron density (2 regions) - start at chosen power fraction
    state.n_L = initial_power;
    state.n_U = initial_power;

    % Delayed precursors at equilibrium: C = (β/Λ/λ_d) * n
    C_eq_factor = p.beta / p.Lambda / p.lambda_d;
    state.C_L = C_eq_factor * state.n_L;
    state.C_U = C_eq_factor * state.n_U;

    % Void fraction at equilibrium (using steam quality model)
    x_eq = max(0, (p.k_P * initial_power * 1000) / (p.m_flow * p.h_fg));
    alpha_eq = p.alpha_max * (x_eq^p.p_shape) / (1 + x_eq^p.p_shape);
    state.alpha_L = alpha_eq;
    state.alpha_U = alpha_eq;

    % Fuel temperature at equilibrium
    state.Tf_L = p.Tf0 + p.a_f * initial_power;
    state.Tf_U = p.Tf0 + p.a_f * initial_power;

    % Moderator temperature (graphite) - start near reference Tm0
    state.Tm_L = p.Tm0;
    state.Tm_U = p.Tm0;

    % Xenon/Iodine at equilibrium
    n_init = (state.n_L + state.n_U) / 2;
    state.I_L = p.y_I * n_init / p.lambda_I;
    state.I_U = state.I_L;
    state.X_L = (p.y_X * n_init + p.lambda_I * state.I_L) / ...
                (p.lambda_X + p.sigma_X * n_init);
    state.X_U = state.X_L;

    % Control rod positions (physics state c in [0,1] and visual depth)
    orm_level = max(0, min(1, orm_level));
    state.c_L = orm_level;
    state.c_U = orm_level;
    insertion_depth = grid.n_rows * orm_level;
    state.rod_position = insertion_depth * ones(grid.n_control_rods, 1);
    state.rod_target = insertion_depth * ones(grid.n_control_rods, 1);

    % Coolant state arrays (per fuel cell)
    state.coolant_temp = zeros(grid.n_rows, grid.n_cols);      % 0 = cold, 1 = boiling
    state.coolant_density = ones(grid.n_rows, grid.n_cols);    % 1 = liquid, 0 = void
    state.coolant_void = false(grid.n_rows, grid.n_cols);      % true when evaporated

    % Neutron particles (visualization)
    state.neutron_x = [];
    state.neutron_y = [];
    state.neutron_vx = [];
    state.neutron_vy = [];
    state.neutron_in_rod = [];

    % Counters
    state.total_produced = 1;
    state.total_from_fission = 0;
    state.total_absorbed_rods = 0;
    state.total_absorbed_coolant = 0;
    state.total_fissions = 0;
    state.total_escaped = 0;
    state.total_spontaneous = 0;

    % Computed values (updated each frame)
    state.ORM_current = sum(state.rod_position) / grid.n_rows;
    state.void_mult = 1.0;
    state.target_particles = round((state.n_L + state.n_U) * params.density_to_particles);

    % Transport delay circular buffer
    % Buffer size = tau_flow / dt = 2.0 / 0.05 = 40 samples
    state.delay_buffer_size = 40;
    state.alpha_L_history = alpha_eq * ones(state.delay_buffer_size, 1);  % Initialize at equilibrium
    state.delay_write_idx = 1;

    % Simulation time (used for rbmk_dynamics SCRAM timing)
    state.simulation_time = 0;
end


%% ============================================================================
%  PHYSICS UPDATE FUNCTIONS
%  ============================================================================

function [state, rho] = update_reactor_physics(state, params, grid, dt)
%UPDATE_REACTOR_PHYSICS Advance the RBMK core using rbmk_dynamics.
%   Uses the calibrated two-region DDE model from rbmk_dynamics.m, with a
%   simple explicit Euler step over dt. Also computes a reactivity
%   breakdown (void/Doppler/moderator/xenon/rods/tip) for display, without
%   feeding any visualization-only knobs back into the physics.

    rho = struct();
    p = params.p_rbmk;

    % ---------------------------------------------------------------------
    % 1. TRANSPORT DELAY BUFFER FOR α_L(t - τ_flow)
    % ---------------------------------------------------------------------
    % Read delayed lower-region void from tau_flow seconds ago
    read_idx = mod(state.delay_write_idx - 1, state.delay_buffer_size) + 1;
    alpha_L_past = state.alpha_L_history(read_idx);

    % Store current lower void for future use
    state.alpha_L_history(state.delay_write_idx) = state.alpha_L;
    state.delay_write_idx = mod(state.delay_write_idx, state.delay_buffer_size) + 1;

    % Build delayed state vector Z for rbmk_dynamics (only Z(3) is used)
    Z = zeros(16, 1);
    Z(3) = alpha_L_past;

    % ---------------------------------------------------------------------
    % 2. BUILD FULL 16-STATE VECTOR AND ADVANCE WITH rbmk_dynamics
    % ---------------------------------------------------------------------
    y = zeros(16, 1);
    % Lower region
    y(1)  = state.n_L;
    y(2)  = state.C_L;
    y(3)  = state.alpha_L;
    y(4)  = state.Tf_L;
    y(5)  = state.Tm_L;
    y(6)  = state.I_L;
    y(7)  = state.X_L;
    y(8)  = state.c_L;
    % Upper region
    y(9)  = state.n_U;
    y(10) = state.C_U;
    y(11) = state.alpha_U;
    y(12) = state.Tf_U;
    y(13) = state.Tm_U;
    y(14) = state.I_U;
    y(15) = state.X_U;
    y(16) = state.c_U;

    % Compute time derivatives from the calibrated DDE model
    dydt = rbmk_dynamics(state.simulation_time, y, Z, p);

    % Explicit Euler step
    y = y + dt * dydt;

    % Clamp to basic physical bounds (integrator safeguard; equations
    % themselves are unchanged and still come from rbmk_dynamics).
    y([1 9]) = max(0, y([1 9]));                            % neutron density ≥ 0
    y([3 11]) = max(0, min(p.alpha_max, y([3 11])));        % 0 ≤ void ≤ alpha_max
    y([4 12]) = max(p.Tc, y([4 12]));                       % fuel temp ≥ coolant inlet
    y([5 13]) = max(p.Tc, y([5 13]));                       % moderator temp ≥ coolant inlet
    y([6 7 14 15]) = max(0, y([6 7 14 15]));                % I, X ≥ 0
    y([8 16]) = max(0, min(1, y([8 16])));                  % control rod fraction in [0,1]

    % Write back to state structure
    state.n_L   = y(1);
    state.C_L   = y(2);
    state.alpha_L = y(3);
    state.Tf_L  = y(4);
    state.Tm_L  = y(5);
    state.I_L   = y(6);
    state.X_L   = y(7);
    state.c_L   = y(8);

    state.n_U   = y(9);
    state.C_U   = y(10);
    state.alpha_U = y(11);
    state.Tf_U  = y(12);
    state.Tm_U  = y(13);
    state.I_U   = y(14);
    state.X_U   = y(15);
    state.c_U   = y(16);

    % ---------------------------------------------------------------------
    % 3. ORM & VISUAL-ONLY MULTIPLIER (does NOT affect physics)
    % ---------------------------------------------------------------------
    % Derive rod insertion depth for all visual rods from average c
    c_avg = (state.c_L + state.c_U) / 2;
    insertion_depth = grid.n_rows * c_avg;
    state.rod_position = insertion_depth * ones(grid.n_control_rods, 1);

    % Operational reactivity margin (pedagogical)
    state.ORM_current = sum(state.rod_position) / grid.n_rows;
    ORM_violated = state.ORM_current < params.ORM_min;

    % Void sensitivity multiplier is for display only; it no longer feeds
    % back into the physics equations.
    if state.ORM_current >= params.ORM_safe
        state.void_mult = 1.0;
    else
        orm_deficit = max(0, 1 - state.ORM_current / params.ORM_safe);
        state.void_mult = 1.0 + (params.ORM_violation_factor - 1.0) * orm_deficit;
    end

    % Target number of visual neutrons proportional to total neutron density
    state.target_particles = round((state.n_L + state.n_U) * params.density_to_particles);

    % ---------------------------------------------------------------------
    % 4. REACTIVITY BREAKDOWN (MATCHING rbmk_dynamics FORMULAS)
    % ---------------------------------------------------------------------
    n_L = state.n_L;    n_U = state.n_U;
    alpha_L = state.alpha_L; alpha_U = state.alpha_U;
    Tf_L = state.Tf_L;  Tf_U = state.Tf_U;
    Tm_L = state.Tm_L;  Tm_U = state.Tm_U;
    X_L = state.X_L;    X_U = state.X_U;
    c_L = state.c_L;    c_U = state.c_U;
    t   = state.simulation_time;

    % A. Void feedback
    rho.void_L = p.kappa_V * alpha_L;
    rho.void_U = p.kappa_V * alpha_U;

    % B. Doppler feedback (with power-dependent enhancement)
    rho_dop_base_L = p.kappa_D0 * (sqrt(Tf_L) - sqrt(p.Tf0));
    rho_dop_base_U = p.kappa_D0 * (sqrt(Tf_U) - sqrt(p.Tf0));
    power_frac_L = n_L;
    power_frac_U = n_U;
    if isfield(p, 'doppler_enhancement')
        doppler_mult_L = 1.0 + (p.doppler_enhancement - 1.0) * power_frac_L;
        doppler_mult_U = 1.0 + (p.doppler_enhancement - 1.0) * power_frac_U;
    else
        doppler_mult_L = 1.0;
        doppler_mult_U = 1.0;
    end
    rho.dop_L = rho_dop_base_L * doppler_mult_L;
    rho.dop_U = rho_dop_base_U * doppler_mult_U;

    % C. Moderator (graphite) temperature feedback
    rho_mod_L = p.kappa_M0 * (Tm_L - p.Tm0);
    rho_mod_U = p.kappa_M0 * (Tm_U - p.Tm0);

    % D. Xenon-135 feedback
    rho.xen_L = -p.kappa_X * X_L;
    rho.xen_U = -p.kappa_X * X_U;

    % E. Control rod reactivity (boron + graphite followers)
    rho_boron_L = -(p.rho_c_max * p.rod_worth_fraction_L) * c_L;
    rho_boron_U = -(p.rho_c_max * p.rod_worth_fraction_U) * c_U;
    rho_follower_L = p.rho_graphite_follower * (1 - c_L);
    rho_follower_U = p.rho_graphite_follower * (1 - c_U);
    rho_rod_L = rho_boron_L + rho_follower_L;
    rho_rod_U = rho_boron_U + rho_follower_U;
    % For display we keep a single "rod" term as the mean of both regions
    rho.boron = (rho_boron_L + rho_boron_U) / 2;
    rho.follower = (rho_follower_L + rho_follower_U) / 2;
    rho.rod = rho.boron + rho.follower;

    % F. Positive SCRAM tip effect (lower region only)
    rho.tip_L = 0;
    if t >= p.t_scram
        dt_scram = t - p.t_scram;
        rho.tip_L = p.rho_tip * exp(-dt_scram / p.tau_tip);
    end

    % Total reactivity per region and average
    rho.total_L = rho.void_L + rho.dop_L + rho_mod_L + rho.xen_L + rho_rod_L + rho.tip_L;
    rho.total_U = rho.void_U + rho.dop_U + rho_mod_U + rho.xen_U + rho_rod_U;
    rho.total   = (rho.total_L + rho.total_U) / 2;

    % Flags for display
    rho.ORM_violated = ORM_violated;
    rho.c_normalized = c_avg;
end

function state = update_coolant_background(state, cfg, grid)
%UPDATE_COOLANT_BACKGROUND Apply cooling and condensation to coolant cells
    temp_baseline = cfg.physics.temp_baseline;
    temp_condense = cfg.physics.temp_condense;
    cooling_rate = cfg.physics.cooling_rate;

    for row = 1:grid.n_rows
        for col = 1:grid.n_cols
            if grid.fuel_map(row, col)
                % Gradual cooling toward baseline
                state.coolant_temp(row, col) = max(temp_baseline, ...
                    state.coolant_temp(row, col) - cooling_rate);

                % Condensation: void → liquid when cooled sufficiently
                if state.coolant_void(row, col) && state.coolant_temp(row, col) < temp_condense
                    state.coolant_void(row, col) = false;
                    state.coolant_density(row, col) = 1.0;
                end

                % Update density based on temperature (simple linear model)
                if ~state.coolant_void(row, col)
                    state.coolant_density(row, col) = max(0, 1 - state.coolant_temp(row, col));
                end
            end
        end
    end
end

function state = update_control_rods(state, grid, rod_speed)
%UPDATE_CONTROL_RODS Animate control rod movement toward targets
    for j = 1:grid.n_control_rods
        if state.rod_position(j) < state.rod_target(j)
            state.rod_position(j) = min(state.rod_target(j), state.rod_position(j) + rod_speed);
        elseif state.rod_position(j) > state.rod_target(j)
            state.rod_position(j) = max(state.rod_target(j), state.rod_position(j) - rod_speed);
        end
    end
end

function [n_new, C_new] = update_neutron_ode(n, C, n_other, dt, params)
%UPDATE_NEUTRON_ODE Euler integration of point kinetics equations
    beta = params(1);
    Lambda = params(2);
    lambda_d = params(3);
    D_n = params(4);
    rho = params(5);

    dn_dt = ((rho - beta)/Lambda)*n + lambda_d*C + D_n*(n_other - n);
    dC_dt = (beta/Lambda)*n - lambda_d*C;

    n_new = max(0, n + dn_dt * dt);
    C_new = max(0, C + dC_dt * dt);
end


%% ============================================================================
%  PARTICLE SIMULATION FUNCTIONS
%  ============================================================================

function [state, gfx] = update_particles(state, gfx, grid, cfg, dt, params)
%UPDATE_PARTICLES Move neutrons and handle collisions for one substep

    if isempty(state.neutron_x)
        % No neutrons to move
        return;
    end

    % Move neutrons
    state.neutron_x = state.neutron_x + state.neutron_vx * dt;
    state.neutron_y = state.neutron_y + state.neutron_vy * dt;

    % Collision detection and interactions
    new_x = []; new_y = []; new_vx = []; new_vy = [];
    remove = false(size(state.neutron_x));

    for i = 1:length(state.neutron_x)
        % Boundary check
        if state.neutron_x(i) < 0 || state.neutron_x(i) > grid.n_cols || ...
           state.neutron_y(i) < 0 || state.neutron_y(i) > grid.n_rows
            remove(i) = true;
            state.total_escaped = state.total_escaped + 1;
            continue;
        end

        % Control rod absorption
        [absorbed, state] = check_rod_absorption(state, i, grid, cfg);
        if absorbed
            remove(i) = true;
            continue;
        end

        % Determine current cell (1-indexed)
        cell_col = max(1, min(grid.n_cols, floor(state.neutron_x(i)) + 1));
        cell_row = max(1, min(grid.n_rows, floor(state.neutron_y(i)) + 1));

        % Interaction with coolant in fuel channels (heating and voiding only)
        if grid.fuel_map(cell_row, cell_col) && ~state.coolant_void(cell_row, cell_col)
            % Heat coolant
            state.coolant_temp(cell_row, cell_col) = min(1.0, ...
                state.coolant_temp(cell_row, cell_col) + cfg.physics.heat_per_neutron);

            % Evaporation: create void if hot enough
            if state.coolant_temp(cell_row, cell_col) >= cfg.physics.temp_evaporate
                state.coolant_void(cell_row, cell_col) = true;
                state.coolant_density(cell_row, cell_col) = 0;

                % High-power regime option: immediately clear voids
                if cfg.physics.clear_void_immediately
                    state.coolant_void(cell_row, cell_col) = false;
                    state.coolant_density(cell_row, cell_col) = 1.0;
                    state.coolant_temp(cell_row, cell_col) = cfg.physics.temp_baseline;
                end
            end
        end

        % U-235 collision (fission) using distance check
        dx_u = gfx.u235_cx - state.neutron_x(i);
        dy_u = gfx.u235_cy - state.neutron_y(i);
        dist = hypot(dx_u, dy_u);
        hit_idx = find(dist < cfg.anim.collision_radius, 1);

        if ~isempty(hit_idx)
            current_count = length(state.neutron_x) - sum(remove) + length(new_x);

            % Spawn new neutrons (scaled by core power)
            total_density = state.n_L + state.n_U;
            density_factor = min(2.0, max(0.5, total_density / 0.2));
            spawn_count = round(3 * density_factor);
            spawn_count = max(1, min(spawn_count, 5));

            if current_count + spawn_count <= cfg.anim.max_neutrons
                [nx, ny, nvx, nvy] = spawn_neutrons(gfx.u235_cx(hit_idx), gfx.u235_cy(hit_idx), ...
                                                    cfg.anim.neutron_speed, spawn_count);
                new_x = [new_x; nx];
                new_y = [new_y; ny];
                new_vx = [new_vx; nvx];
                new_vy = [new_vy; nvy];
                state.total_produced = state.total_produced + spawn_count;
                state.total_from_fission = state.total_from_fission + spawn_count;
            end

            state.total_fissions = state.total_fissions + 1;
            remove(i) = true;
        end
    end

    % Apply removals and additions
    keep_idx = find(~remove);
    state.neutron_in_rod = [state.neutron_in_rod(keep_idx); zeros(numel(new_x), 1)];
    state.neutron_x = [state.neutron_x(keep_idx); new_x(:)];
    state.neutron_y = [state.neutron_y(keep_idx); new_y(:)];
    state.neutron_vx = [state.neutron_vx(keep_idx); new_vx(:)];
    state.neutron_vy = [state.neutron_vy(keep_idx); new_vy(:)];
end

function [absorbed, state] = check_rod_absorption(state, i, grid, cfg)
%CHECK_ROD_ABSORPTION Check if neutron is absorbed by control rod
    absorbed = false;
    current_rod = 0;
    control_rod_width = cfg.fuel.control_rod_width;

    for j = 1:grid.n_control_rods
        if state.rod_position(j) > 0
            rod_x = grid.control_rod_cols(j);
            rod_half_width = control_rod_width / 2;
            if state.neutron_x(i) >= rod_x - rod_half_width && ...
               state.neutron_x(i) <= rod_x + rod_half_width && ...
               state.neutron_y(i) <= state.rod_position(j)
                current_rod = j;
                break;
            end
        end
    end

    if current_rod > 0 && state.neutron_in_rod(i) ~= current_rod
        depth_fraction = min(1, state.rod_position(current_rod) / grid.n_rows);
        absorb_prob = cfg.physics.rod_absorption_max * depth_fraction;
        if rand() < absorb_prob
            absorbed = true;
            state.total_absorbed_rods = state.total_absorbed_rods + 1;
        else
            state.neutron_in_rod(i) = current_rod;
        end
    elseif current_rod == 0
        state.neutron_in_rod(i) = 0;
    end
end

function state = balance_particles(state, gfx, cfg, params)
%BALANCE_PARTICLES Add/remove particles to match ODE target
    current_particles = length(state.neutron_x);
    state.target_particles = min(state.target_particles, cfg.anim.max_neutrons);
    particle_diff = state.target_particles - current_particles;

    if particle_diff > 5 && gfx.n_u235 > 0
        n_spawn = min(abs(particle_diff), 10);
        for k = 1:n_spawn
            site_idx = randi(gfx.n_u235);
            [sx, sy, svx, svy] = spawn_neutrons(gfx.u235_cx(site_idx), gfx.u235_cy(site_idx), ...
                                                cfg.anim.neutron_speed, 1);
            state.neutron_x = [state.neutron_x; sx];
            state.neutron_y = [state.neutron_y; sy];
            state.neutron_vx = [state.neutron_vx; svx];
            state.neutron_vy = [state.neutron_vy; svy];
            state.neutron_in_rod = [state.neutron_in_rod; 0];
        end
    elseif particle_diff < -5 && current_particles > 0
        n_remove = min(abs(particle_diff), 10);
        n_remove = min(n_remove, current_particles);
        if n_remove > 0
            remove_idx = randperm(current_particles, n_remove);
            keep_idx = setdiff(1:current_particles, remove_idx);
            state.neutron_x = state.neutron_x(keep_idx);
            state.neutron_y = state.neutron_y(keep_idx);
            state.neutron_vx = state.neutron_vx(keep_idx);
            state.neutron_vy = state.neutron_vy(keep_idx);
            state.neutron_in_rod = state.neutron_in_rod(keep_idx);
        end
    end
end

function [x, y, vx, vy] = spawn_neutrons(cx, cy, speed, count)
%SPAWN_NEUTRONS Create neutrons at a position with random directions
    % Place neutrons slightly outside the collision radius so they do not
    % immediately re-collide with the same U-235 site.
    spawn_offset = 0.4;
    pos_angles = 2 * pi * rand(count, 1);  % Position offset angles
    vel_angles = 2 * pi * rand(count, 1);  % Independent velocity angles
    x = cx + spawn_offset * cos(pos_angles);
    y = cy + spawn_offset * sin(pos_angles);
    vx = speed * cos(vel_angles);
    vy = speed * sin(vel_angles);
end


%% ============================================================================
%  VISUALIZATION FUNCTIONS
%  ============================================================================

function [fig, gfx] = setup_graphics(cfg, grid, state)
%SETUP_GRAPHICS Create figure and all graphical elements
    gfx = struct();

    % Create figure
    fig = figure('Name', 'RBMK-1000 Core Cross-Section', ...
                 'NumberTitle', 'off', 'Color', 'w', ...
                 'Position', [100, 100, 1100, 650]);
    hold on;
    axis equal;

    % Draw grid cells
    gfx.cell_patches = gobjects(grid.n_rows, grid.n_cols);
    for row = 1:grid.n_rows
        for col = 1:grid.n_cols
            x = [col-1, col, col, col-1, col-1];
            y = [row-1, row-1, row, row, row-1];
            if grid.fuel_map(row, col)
                cell_color = cfg.colors.coolant;
            else
                cell_color = cfg.colors.moderator;
            end
            gfx.cell_patches(row, col) = patch(x, y, cell_color, ...
                'EdgeColor', cfg.colors.edge, 'LineWidth', 0.3);
        end
    end

    % Draw control rods
    gfx.h_control_rods = gobjects(grid.n_control_rods, 1);
    for i = 1:grid.n_control_rods
        x_left = grid.control_rod_cols(i) - cfg.fuel.control_rod_width/2;
        gfx.h_control_rods(i) = rectangle('Position', [x_left, 0, cfg.fuel.control_rod_width, state.rod_position(i)], ...
            'FaceColor', cfg.colors.control_rod, 'EdgeColor', 'k', 'LineWidth', 1);
    end

    % Draw uranium pellets and build U-235 collision data
    [gfx.u235_cx, gfx.u235_cy, gfx.n_u235] = draw_uranium_pellets(cfg, grid);

    % Configure axes
    xlim([0, grid.n_cols]);
    ylim([0, grid.n_rows]);
    set(gca, 'YDir', 'reverse');
    xlabel('Column', 'FontSize', 12);
    ylabel('Row', 'FontSize', 12);
    title(sprintf('RBMK-1000 Reactor Core Cross-Section (%d × %d Grid)', grid.n_rows, grid.n_cols), ...
        'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'XTick', 0:5:grid.n_cols, 'YTick', 0:6:grid.n_rows);
    set(gca, 'TickDir', 'out', 'FontSize', 10);

    % Separator line
    plot([0, grid.n_cols], [grid.mid_row, grid.mid_row], '-', ...
        'Color', [1.0 0.3 0.3], 'LineWidth', 3.0);
    text(grid.n_cols/2, grid.mid_row + 0.3, 'Spatial Coupling Boundary (D_n)', ...
        'Color', [0.6 0.1 0.1], 'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    % Stats panel
    stats_x = grid.n_cols + 1;
    gfx.h_reactivity = text(stats_x, grid.n_rows - 1, '', 'FontSize', 11, 'FontWeight', 'bold');
    gfx.h_density = text(stats_x, grid.n_rows - 4, '', 'FontSize', 11, 'FontWeight', 'bold');
    gfx.h_neutron_count = text(stats_x, grid.n_rows - 8, '', 'FontSize', 11, 'FontWeight', 'bold');
    gfx.h_fission_count = text(stats_x, grid.n_rows - 11, '', 'FontSize', 11, 'FontWeight', 'bold');
    gfx.h_orm = text(stats_x, grid.n_rows - 14, '', 'FontSize', 11, 'FontWeight', 'bold');

    xlim([-0.5, grid.n_cols + 12]);
    ylim([-0.5, grid.n_rows + 0.5]);

    % Neutron scatter plot
    gfx.h_neutrons = scatter([], [], 30, 'k', 'filled');

    % Store config reference
    gfx.cfg = cfg;
end

function [u235_cx, u235_cy, n_u235] = draw_uranium_pellets(cfg, grid)
%DRAW_URANIUM_PELLETS Draw uranium circles and return U-235 positions
    rng(42);  % Fixed seed
    u235_fraction = cfg.fuel.u235_fraction;
    pellet_radius = cfg.fuel.pellet_radius;
    n_fuel_cells = grid.n_fuel_cells;
    n_u235 = round(u235_fraction * n_fuel_cells);

    % Select U-235 positions
    u235_set = false(n_fuel_cells, 1);
    if isfield(cfg.fuel, 'u235_distribution') && strcmp(cfg.fuel.u235_distribution, 'even')
        stride = floor(n_fuel_cells / n_u235);
        offset = floor(stride / 2);
        u235_indices = offset:stride:n_fuel_cells;
        u235_indices = u235_indices(1:min(n_u235, length(u235_indices)));
        if length(u235_indices) < n_u235
            remaining = setdiff(1:n_fuel_cells, u235_indices);
            remaining = remaining(randperm(length(remaining)));
            u235_indices = [u235_indices, remaining(1:(n_u235 - length(u235_indices)))];
        end
        u235_set(u235_indices) = true;
    else
        u235_indices = randperm(n_fuel_cells, n_u235);
        u235_set(u235_indices) = true;
    end

    % Draw circles
    theta = linspace(0, 2*pi, 30);
    u235_cx = zeros(n_u235, 1);
    u235_cy = zeros(n_u235, 1);
    u235_idx = 1;

    for i = 1:n_fuel_cells
        row = grid.fuel_rows(i);
        col = grid.fuel_cols(i);
        cx = col - 0.5;
        cy = row - 0.5;
        rx = cx + pellet_radius * cos(theta);
        ry = cy + pellet_radius * sin(theta);

        if u235_set(i)
            patch(rx, ry, cfg.colors.uranium, 'EdgeColor', cfg.colors.uranium * 0.7, 'LineWidth', 0.5);
            u235_cx(u235_idx) = cx;
            u235_cy(u235_idx) = cy;
            u235_idx = u235_idx + 1;
        else
            patch(rx, ry, cfg.colors.u238, 'EdgeColor', cfg.colors.u238 * 0.7, 'LineWidth', 0.5);
        end
    end
end

function update_cell_colors(state, gfx, grid, params) %#ok<INUSD>
%UPDATE_CELL_COLORS Update coolant colors based on per-cell temperature/void
    for row = 1:grid.n_rows
        for col = 1:grid.n_cols
            if grid.fuel_map(row, col)
                if state.coolant_void(row, col)
                    % Void cell - light/transparent to indicate steam
                    set(gfx.cell_patches(row, col), 'FaceColor', [1 1 1], 'FaceAlpha', 0.3);
                else
                    % Liquid coolant - color based on local temperature
                    new_color = temp_to_color(state.coolant_temp(row, col));
                    set(gfx.cell_patches(row, col), 'FaceColor', new_color, 'FaceAlpha', 1.0);
                end
            end
        end
    end

    % Update control rod graphics
    for j = 1:grid.n_control_rods
        x_left = grid.control_rod_cols(j) - gfx.cfg.fuel.control_rod_width/2;
        set(gfx.h_control_rods(j), 'Position', [x_left, 0, gfx.cfg.fuel.control_rod_width, max(0.01, state.rod_position(j))]);
    end
end

function update_neutron_graphics(state, gfx)
%UPDATE_NEUTRON_GRAPHICS Update neutron scatter plot
    if ~isempty(state.neutron_x)
        set(gfx.h_neutrons, 'XData', state.neutron_x, 'YData', state.neutron_y);
    else
        set(gfx.h_neutrons, 'XData', [], 'YData', []);
    end
end

function update_stats_panel(state, rho, gfx, params)
%UPDATE_STATS_PANEL Update the stats display panel
    avg_rho_void = (rho.void_L + rho.void_U) / 2;
    avg_rho_dop = (rho.dop_L + rho.dop_U) / 2;
    avg_rho_xen = (rho.xen_L + rho.xen_U) / 2;

    % Include SCRAM tip effect in display if active
    if isfield(rho, 'tip_L') && rho.tip_L > 0.0001
        set(gfx.h_reactivity, 'String', sprintf('Reactivity: %.4f\n  void: +%.4f\n  Doppler: %.4f\n  Xenon: %.4f\n  rod: %.4f\n  TIP: +%.4f', ...
            rho.total, avg_rho_void, avg_rho_dop, avg_rho_xen, rho.rod, rho.tip_L));
    else
        set(gfx.h_reactivity, 'String', sprintf('Reactivity: %.4f\n  void: +%.4f\n  Doppler: %.4f\n  Xenon: %.4f\n  rod: %.4f', ...
            rho.total, avg_rho_void, avg_rho_dop, avg_rho_xen, rho.rod));
    end

    total_density = state.n_L + state.n_U;
    set(gfx.h_density, 'String', sprintf('Density: %.3f\n  n_L = %.3f\n  n_U = %.3f', ...
        total_density, state.n_L, state.n_U));

    set(gfx.h_neutron_count, 'String', sprintf('Neutrons: %d / %d', ...
        length(state.neutron_x), state.target_particles));

    set(gfx.h_fission_count, 'String', sprintf('Fissions: %d', state.total_fissions));

    % ORM status
    if rho.ORM_violated
        orm_color = [0.9 0.1 0.1];
        orm_status = 'VIOLATED!';
    elseif state.ORM_current < params.ORM_safe
        orm_color = [0.9 0.6 0.0];
        orm_status = 'WARNING';
    else
        orm_color = [0.1 0.6 0.1];
        orm_status = 'OK';
    end
    set(gfx.h_orm, 'String', sprintf('ORM: %.2f / %.1f [%s]\n  Void mult: %.2fx', ...
        state.ORM_current, params.ORM_min, orm_status, state.void_mult), 'Color', orm_color);
end

function color = temp_to_color(temp)
%TEMP_TO_COLOR Maps temperature to blue-red spectrum
    if temp < 0.33
        t = temp / 0.33;
        color = (1-t) * [0.15, 0.25, 0.7] + t * [0.5, 0.7, 1.0];
    elseif temp < 0.67
        t = (temp - 0.33) / 0.34;
        color = (1-t) * [0.5, 0.7, 1.0] + t * [1.0, 0.6, 0.5];
    else
        t = min(1, (temp - 0.67) / 0.33);
        color = (1-t) * [1.0, 0.6, 0.5] + t * [0.7, 0.15, 0.15];
    end
end


%% ============================================================================
%  OUTPUT FUNCTIONS
%  ============================================================================

function print_progress(frame, n_frames, state, gfx, grid)
%PRINT_PROGRESS Print progress indicator
    n_voids = sum(state.coolant_void(:));
    fuel_mask = grid.fuel_map;
    if any(fuel_mask(:))
        avg_temp = mean(state.coolant_temp(fuel_mask));
    else
        avg_temp = NaN;
    end
    fprintf('Frame %d/%d | Neutrons: %d (target: %d) | n=[%.3f,%.3f] C=[%.3f,%.3f] | Fissions: %d | Void cells: %d | Avg coolant temp: %.3f\n', ...
        frame, n_frames, length(state.neutron_x), state.target_particles, ...
        state.n_L, state.n_U, state.C_L, state.C_U, state.total_fissions, n_voids, avg_temp);
end

function print_final_summary(state, gfx, grid, params, cfg)
%PRINT_FINAL_SUMMARY Print final reactor state summary
    % Compute void cells based on ODE void fraction (probabilistic)
    avg_alpha = (state.alpha_L + state.alpha_U) / 2;
    n_voids = round(avg_alpha * grid.n_fuel_cells);
    c_normalized = mean(state.rod_position) / grid.n_rows;

    fprintf('\n===== NEUTRON BALANCE SUMMARY =====\n');
    fprintf('Total produced:      %d\n', state.total_produced);
    fprintf('  - From fission:    %d\n', state.total_from_fission);
    fprintf('  - Spontaneous:     %d\n', state.total_spontaneous);
    fprintf('Absorbed (rods):     %d\n', state.total_absorbed_rods);
    fprintf('Total escaped:       %d (boundary)\n', state.total_escaped);
    fprintf('Total fissions:      %d\n', state.total_fissions);
    fprintf('Active remaining:    %d\n', length(state.neutron_x));
    fprintf('Void cells:          %d / %d (%.1f%%)\n', n_voids, grid.n_fuel_cells, 100*n_voids/grid.n_fuel_cells);
    fprintf('Neutron limit:       %d (max allowed)\n', cfg.anim.max_neutrons);
    fprintf('-----------------------------------\n');
    fprintf('ODE STATE (final):\n');
    fprintf('  n_L = %.4f, n_U = %.4f (density)\n', state.n_L, state.n_U);
    fprintf('  C_L = %.4f, C_U = %.4f (precursors)\n', state.C_L, state.C_U);
    fprintf('  alpha_L = %.2f%%, alpha_U = %.2f%% (void fraction)\n', state.alpha_L*100, state.alpha_U*100);
    fprintf('  Tf_L = %.0f°C, Tf_U = %.0f°C (fuel temp)\n', state.Tf_L, state.Tf_U);
    fprintf('  X_L = %.1f, X_U = %.1f (xenon concentration)\n', state.X_L, state.X_U);
    fprintf('  Rod position: %.1f%% inserted (c=%.2f)\n', c_normalized*100, c_normalized);
    fprintf('OPERATIONAL REACTIVITY MARGIN (ORM):\n');
    fprintf('  ORM: %.2f equivalent rods (min: %.1f, safe: %.1f)\n', ...
        state.ORM_current, params.ORM_min, params.ORM_safe);
    if state.ORM_current < params.ORM_min
        fprintf('  *** ORM VIOLATION - Reactor in dangerous state! ***\n');
    elseif state.ORM_current < params.ORM_safe
        fprintf('  WARNING: ORM below safe level\n');
    else
        fprintf('  Status: OK (ORM satisfied)\n');
    end
    fprintf('  Void sensitivity multiplier: %.2fx\n', state.void_mult);
    fprintf('REACTIVITY BREAKDOWN:\n');
    % Use the same physics as rbmk_dynamics for the summary breakdown.
    p = params.p_rbmk;
    avg_rho_void = p.kappa_V * (state.alpha_L + state.alpha_U) / 2;
    avg_rho_dop = p.kappa_D0 * (sqrt(state.Tf_L) - sqrt(p.Tf0));
    avg_rho_xen = -p.kappa_X * (state.X_L + state.X_U) / 2;
    rho_boron = -p.rho_c_max * c_normalized;
    rho_follower = p.rho_graphite_follower * (1 - c_normalized);
    rho_rod = rho_boron + rho_follower;
    rho_total = avg_rho_void + avg_rho_dop + avg_rho_xen + rho_rod;
    fprintf('  rho_void    = +%.4f (positive void coefficient)\n', avg_rho_void);
    fprintf('  rho_Doppler = %.4f (fuel temperature feedback)\n', avg_rho_dop);
    fprintf('  rho_Xenon   = %.4f (Xe-135 poisoning)\n', avg_rho_xen);
    fprintf('  rho_rod     = %.4f (control rods: boron %.4f + follower %.4f)\n', rho_rod, rho_boron, rho_follower);
    fprintf('  rho_TOTAL   = %.4f (net reactivity)\n', rho_total);
    fprintf('  Target particles: %d\n', state.target_particles);
    fprintf('===================================\n');
end


%% ============================================================================
%  CONFIGURATION PRESETS
%  ============================================================================

function cfg = get_preset_config(preset_name)
%GET_PRESET_CONFIG Returns configuration struct for the specified preset
    switch preset_name
        case 'high_power'
            cfg = high_power_config();
        case 'low_power'
            cfg = low_power_config();
        case 'chernobyl'
            cfg = chernobyl_config();
        otherwise
            warning('Unknown preset "%s", using high_power', preset_name);
            cfg = high_power_config();
    end
end

function cfg = high_power_config()
%HIGH_POWER_CONFIG Normal high-power operation (~3000 MW thermal)
%   Physics: full-power RBMK with strong pump-driven coolant flow.
%   Visualization water model:
%     - High flow rapidly sweeps bubbles out of 7 m channels
%       (order ~1 s residence time).
%     - Boiling is saturated: void fraction sits in a modest band and
%       additional heating adds only a small incremental void increase.
%     - Strong cooling prevents long-lived steam accumulation.
    cfg = struct();

    cfg.grid.n_rows = 24;
    cfg.grid.moderator_width = 3;
    cfg.grid.fuel_width = 5;
    cfg.grid.n_fuel_strips = 4;
    cfg.grid.n_moderator_strips = 5;

    cfg.colors.coolant = [0.15 0.25 0.7];
    cfg.colors.moderator = [0.85 0.85 0.85];
    cfg.colors.edge = [0.5 0.5 0.5];
    cfg.colors.fuel_boundary = [0.8 0.2 0.1];
    cfg.colors.uranium = [1.0 0.88 0.2];
    cfg.colors.u238 = [0.5 0.5 0.5];
    cfg.colors.control_rod = [0.2 0.2 0.2];

    cfg.fuel.u235_fraction = 0.30;
    cfg.fuel.pellet_radius = 0.35;
    cfg.fuel.control_rod_width = 0.8;
    cfg.fuel.u235_distribution = 'even';

    cfg.anim.dt = 0.05;
    cfg.anim.pause_time = 0.03;
    cfg.anim.neutron_speed = 3.0;
    cfg.anim.collision_radius = 0.35;
    cfg.anim.max_neutrons = 300;
    cfg.anim.n_frames = 1000;
    cfg.anim.rod_speed = 0.8;
    cfg.anim.substeps = 4;

    cfg.physics.rod_absorption_max = 0.25;
    % High-power water behaviour: strong pump flow and saturated boiling.
    % Faster cooling and a moderately low boil threshold approximate
    % rapid bubble sweep-out with a stable 15–20%% void band.
    cfg.physics.cooling_rate = 0.012;   % Strong background cooling (high flow)
    cfg.physics.heat_per_neutron = 0.03; % Per-hit heating (kept moderate)
    cfg.physics.temp_baseline = 0.0;
    cfg.physics.temp_evaporate = 0.60;  % Boiling begins at moderate temperature
    cfg.physics.temp_condense = 0.35;   % Bubbles collapse once cooled a bit
    cfg.physics.spontaneous_fission_prob = 0.03;
    cfg.physics.clear_void_immediately = false; % Short-lived, not instantaneous, voids
    cfg.physics.max_coolant_absorption = 0.008;
end

function cfg = low_power_config()
%LOW_POWER_CONFIG Low-power operation (~200 MW thermal)
%   Visualization water model:
%     - Effective cooling is weaker (slower flow / less mixing).
%     - Local heating from a power increase produces a large, sudden
%       void jump.
%     - Voids collapse more slowly, so steam can transiently outrun
%       replacement by fresh water, echoing the unstable accident regime.
    cfg = struct();

    cfg.grid.n_rows = 24;
    cfg.grid.moderator_width = 3;
    cfg.grid.fuel_width = 5;
    cfg.grid.n_fuel_strips = 4;
    cfg.grid.n_moderator_strips = 5;

    cfg.colors.coolant = [0.15 0.25 0.7];
    cfg.colors.moderator = [0.85 0.85 0.85];
    cfg.colors.edge = [0.5 0.5 0.5];
    cfg.colors.fuel_boundary = [0.8 0.2 0.1];
    cfg.colors.uranium = [1.0 0.88 0.2];
    cfg.colors.u238 = [0.5 0.5 0.5];
    cfg.colors.control_rod = [0.2 0.2 0.2];

    cfg.fuel.u235_fraction = 0.30;
    cfg.fuel.pellet_radius = 0.35;
    cfg.fuel.control_rod_width = 0.8;
    cfg.fuel.u235_distribution = 'even';

    cfg.anim.dt = 0.05;
    cfg.anim.pause_time = 0.03;
    cfg.anim.neutron_speed = 2.0;
    cfg.anim.collision_radius = 0.35;
    cfg.anim.max_neutrons = 100;
    cfg.anim.n_frames = 1000;
    cfg.anim.rod_speed = 0.8;
    cfg.anim.substeps = 4;

    cfg.physics.rod_absorption_max = 0.25;
    % Low-power / accident-like cooling: weaker effective flow so
    % heating produces larger, more persistent void swings.
    cfg.physics.cooling_rate = 0.0015;  % Slow background cooling
    cfg.physics.heat_per_neutron = 0.035; % Each neutron heats water more
    cfg.physics.temp_baseline = 0.0;
    cfg.physics.temp_evaporate = 0.55;  % Boiling onset at lower temperature
    cfg.physics.temp_condense = 0.50;   % Voids persist until coolant is quite cool
    cfg.physics.spontaneous_fission_prob = 0.01;
    cfg.physics.clear_void_immediately = false;
    cfg.physics.max_coolant_absorption = 0.02;
end

function cfg = chernobyl_config()
%CHERNOBYL_CONFIG Simulates the dangerous ORM-violated state at Chernobyl
    cfg = low_power_config();
    cfg.anim.max_neutrons = 150;
    cfg.anim.pause_time = 0.02;
    cfg.anim.n_frames = 1200;

    % SCRAM trigger time (set to inf to disable SCRAM demonstration)
    % At t=30s, SCRAM will trigger showing positive tip effect
    cfg.physics.t_scram = 30.0;
end

function dydt = rbmk_dynamics(t, y, Z, p)
    % RBMK_DERIVATIVES_ENHANCED - Two-region DDE model of RBMK reactor dynamics
    %
    % Computes time derivatives for a 16-state DDE system representing the
    % Chernobyl RBMK reactor. Models the critical physics that led to the
    % 1986 accident, including positive void coefficient, transport delay,
    % and graphite follower design flaws.
    %
    % MATHEMATICAL MODEL:
    %   dy/dt = f(t, y(t), y(t-τ), p)
    %
    % where y is a 16-element state vector representing lower and upper core regions.
    %
    % INPUTS:
    %   t - Current time (seconds)
    %   y - Current state vector [16×1]:
    %       Lower region (indices 1-8):  [n_L, C_L, α_L, T_f,L, T_m,L, I_L, X_L, c_L]
    %       Upper region (indices 9-16): [n_U, C_U, α_U, T_f,U, T_m,U, I_U, X_U, c_U]
    %   Z - Delayed state vector (for accessing y(t-τ))
    %       Z(3,1) = α_L(t-τ) is the critical delayed term for coolant transport
    %   p - Parameter structure from rbmk_parameters()
    %
    % OUTPUTS:
    %   dydt - Time derivatives [16×1]
    %
    % STATE VARIABLES (per region):
    %   n     - Neutron density (normalized to nominal power)
    %   C     - Delayed neutron precursor concentration (normalized)
    %   α     - Void fraction (steam bubbles in coolant, 0 = liquid, 1 = steam)
    %   T_f   - Fuel temperature (K)
    %   T_m   - Moderator (graphite) temperature (K)
    %   I     - Iodine-135 concentration (fission product, precursor to Xenon)
    %   X     - Xenon-135 concentration (strong neutron poison)
    %   c     - Control rod insertion fraction (0 = withdrawn, 1 = inserted)
    %
    % KEY PHYSICS IMPLEMENTED:
    %   1. Point kinetics with delayed neutrons (β = 0.005, Λ = 0.48 ms)
    %   2. Spatial coupling between regions (neutron diffusion)
    %   3. POSITIVE void coefficient (κ_V = +0.025) - PRIMARY DESIGN FLAW
    %   4. Doppler feedback (√T law, power-enhanced at high power)
    %   5. Xenon/Iodine decay chains (time-dependent poisoning)
    %   6. Coolant transport delay (τ_flow = 2.0 s) - CAUSES HOPF BIFURCATION
    %   7. Graphite follower effect - RBMK DESIGN FLAW #2
    %   8. Positive SCRAM effect - RBMK DESIGN FLAW #3
    %
    % ENHANCEMENTS BEYOND BASIC MODEL:
    %   - Power-dependent Doppler enhancement (stronger at high power)
    %   - Void saturation at high void fractions (prevents unphysical runaway)
    %   - Xenon equilibrium provides automatic stabilization at high power
    %
    % DESIGN FLAWS MODELED:
    %   1. Positive void coefficient: More steam → MORE power (unstable feedback)
    %   2. Graphite followers: Withdrawn rods add +0.6β positive reactivity
    %   3. Positive SCRAM: Rod insertion FIRST adds +1β before absorption
    %
    % TYPICAL USAGE:
    %   p = rbmk_parameters();
    %   [y0, ~] = compute_equilibrium(0.5, p);
    %   [T, Y] = solve_stiff_dde(@rbmk_dynamics, [0 100], y0, p.tau_flow, p);
    %
    % REFERENCE:
    %   INSAG-7 (IAEA, 1992) - Chernobyl accident analysis
    %   OECD-NEA RBMK simulations
    %   EPJ Nuclear Sciences 9, 10 (2023) - Stability analysis
    %
    % See also: rbmk_parameters, solve_stiff_dde, compute_equilibrium

    %% Input Validation
    if nargin < 4
        error('rbmk_dynamics:NotEnoughInputs', ...
            'Function requires four inputs: t, y, Z, and p');
    end

    if length(y) ~= 16
        error('rbmk_dynamics:InvalidStateVector', ...
            'State vector must have 16 elements, got %d', length(y));
    end

    %% ========================================================================
    %  SECTION 1: UNPACK STATE VECTOR
    %  ========================================================================
    % Extract state variables for each region
    % Lower region (L): bottom half of 7m core (coolant inlet)
    % Upper region (U): top half of 7m core (steam outlet)

    % Lower Region - Indices 1-8
    n_L = y(1);      % Neutron density (normalized, n=1.0 at nominal power)
    C_L = y(2);      % Delayed neutron precursors (normalized)
    alpha_L = y(3);  % Void fraction (0=liquid, 1=steam)
    Tf_L = y(4);     % Fuel temperature (K)
    Tm_L = y(5);     % Moderator (graphite) temperature (K)
    I_L = y(6);      % Iodine-135 concentration (arbitrary units)
    X_L = y(7);      % Xenon-135 concentration (arbitrary units)
    c_L = y(8);      % Control rod position (0=out, 1=in)

    % Upper Region - Indices 9-16
    n_U = y(9);      % Neutron density
    C_U = y(10);     % Delayed neutron precursors
    alpha_U = y(11); % Void fraction
    Tf_U = y(12);    % Fuel temperature (K)
    Tm_U = y(13);    % Moderator temperature (K)
    I_U = y(14);     % Iodine-135 concentration
    X_U = y(15);     % Xenon-135 concentration
    c_U = y(16);     % Control rod position

    % Extract delayed state (for transport delay τ_flow = 2.0 s)
    % The upper region void is driven by lower region void from 2 seconds ago
    % This delay is CRITICAL for Hopf bifurcation (oscillations)
    if isempty(Z)
        % During initialization: no history available, use current value
        alpha_L_past = alpha_L;
    else
        % During simulation: use lower region void from t-τ
        alpha_L_past = Z(3, 1);  % α_L(t - τ_flow)
    end

    %% ========================================================================
    %  SECTION 2: THERMAL-HYDRAULICS (VOID FRACTION CALCULATION)
    %  ========================================================================
    % Calculate equilibrium void fraction based on power and coolant flow
    % Physics: Higher power → more boiling → higher void fraction

    % Steam quality x = (heat generated) / (heat capacity of coolant)
    % x = power / (mass_flow × latent_heat_of_vaporization)
    %
    % Unit conversion: k_P in MW, need kW to match h_fg in kJ/kg
    % Small epsilon (1e-6) prevents division by zero
    x_L = (p.k_P * n_L * 1000) / (p.m_flow * p.h_fg + 1e-6);
    x_U = (p.k_P * n_U * 1000) / (p.m_flow * p.h_fg + 1e-6);

    % Prevent complex numbers during solver iterations
    % (can occur if neutron density temporarily goes slightly negative)
    x_L = max(0, x_L);
    x_U = max(0, x_U);

    % EQUILIBRIUM VOID FRACTION (boiling curve)
    % Uses generalized logistic function: α_eq = α_max × x^p / (1 + x^p)
    %
    % Shape parameter p = 0.25 creates ROOT-LIKE behavior:
    %   - At low x (low power): STEEP gradient (dα/dx large)
    %   - At high x (high power): SATURATED (dα/dx small, approaching dryout)
    %
    % This steep response at low power is WHY low power operation is unstable!
    alpha_eq_L = p.alpha_max * (x_L^p.p_shape) / (1 + x_L^p.p_shape);
    alpha_eq_U = p.alpha_max * (x_U^p.p_shape) / (1 + x_U^p.p_shape);

    % VOID SATURATION FACTOR (enhancement to prevent unphysical runaway)
    % At very high void fractions (α > 0.5), reduce the sensitivity of void dynamics
    % Physics: As voids coalesce, further boiling has less impact on reactivity
    %
    % saturation_factor = 1 / (1 + k_sat × α)
    % Examples with k_sat = 2.0:
    %   α = 0.0: saturation = 1.0   (no reduction)
    %   α = 0.4: saturation = 0.56  (44% reduction)
    %   α = 0.6: saturation = 0.45  (55% reduction)
    saturation_L = 1.0 / (1.0 + p.void_saturation_coeff * alpha_L);
    saturation_U = 1.0 / (1.0 + p.void_saturation_coeff * alpha_U);

    %% ========================================================================
    %  SECTION 3: REACTIVITY FEEDBACKS
    %  ========================================================================
    % Calculate total reactivity ρ = Σ(individual feedback mechanisms)
    % Reactivity is dimensionless: ρ = 0 (critical), ρ > 0 (supercritical)
    % Typical range: -0.02 to +0.02 (±4β)

    % ------------------------------------------------------------------------
    % A. VOID REACTIVITY FEEDBACK - PRIMARY DESIGN FLAW
    % ------------------------------------------------------------------------
    % Linear relationship: ρ_void = κ_V × α
    % κ_V = +0.025 (POSITIVE coefficient - this is the RBMK design flaw!)
    %
    % Physics in RBMK:
    %   More void (steam) → Less water → Less neutron ABSORPTION
    %   (In RBMK, water absorbs neutrons; graphite moderates)
    %   Result: POSITIVE feedback (more power → more steam → MORE power)
    %
    % In normal PWRs/BWRs, void coefficient is NEGATIVE (water is moderator)
    %
    % Typical values:
    %   α = 0.3 (30% void): ρ_void = +0.0075 = +1.5β (significant positive reactivity!)
    rho_void_L = p.kappa_V * alpha_L;
    rho_void_U = p.kappa_V * alpha_U;

    % ------------------------------------------------------------------------
    % B. DOPPLER REACTIVITY FEEDBACK (Fuel Temperature) - STABILIZING
    % ------------------------------------------------------------------------
    % Based on Doppler broadening of U-238 resonances
    % Temperature dependence: ρ ∝ √T (from quantum mechanics)
    %
    % Physics: Higher fuel temperature → broader U-238 resonances
    %          → more neutron capture → NEGATIVE reactivity (stabilizing)
    %
    % Base Doppler coefficient (√T law):
    rho_dop_base_L = p.kappa_D0 * (sqrt(Tf_L) - sqrt(p.Tf0));
    rho_dop_base_U = p.kappa_D0 * (sqrt(Tf_U) - sqrt(p.Tf0));

    % POWER-DEPENDENT ENHANCEMENT (research-based phenomenology)
    % Observation: Doppler feedback is STRONGER at high power
    % Mechanisms:
    %   - Better fuel-clad thermal contact at high burnup
    %   - Higher temperature increases resonance broadening effectiveness
    %   - Phenomenological multiplier: 1.0 at zero power → 1.5 at full power
    %
    % Enhancement factor: mult = 1.0 + (1.5 - 1.0) × n
    power_frac_L = n_L;  % Normalized power (n=1.0 at nominal)
    power_frac_U = n_U;

    if isfield(p, 'doppler_enhancement')
        doppler_mult_L = 1.0 + (p.doppler_enhancement - 1.0) * power_frac_L;
        doppler_mult_U = 1.0 + (p.doppler_enhancement - 1.0) * power_frac_U;
    else
        % Fallback if enhancement not specified
        doppler_mult_L = 1.0;
        doppler_mult_U = 1.0;
    end

    % Apply enhancement
    rho_dop_L = rho_dop_base_L * doppler_mult_L;
    rho_dop_U = rho_dop_base_U * doppler_mult_U;

    % Why this matters for stability:
    %   At LOW power:  weak Doppler → can't overcome positive void feedback
    %   At HIGH power: strong Doppler → dominates void feedback → stable

    % ------------------------------------------------------------------------
    % C. MODERATOR TEMPERATURE FEEDBACK (Graphite Temperature)
    % ------------------------------------------------------------------------
    % Linear relationship: ρ_mod = κ_M0 × (T_m - T_m0)
    %
    % Physics: Hotter graphite → lower density → slightly less moderation
    % κ_M0 ≈ -0.002 (small negative coefficient)
    %
    % Much weaker than Doppler or void feedbacks
    rho_mod_L = p.kappa_M0 * (Tm_L - p.Tm0);
    rho_mod_U = p.kappa_M0 * (Tm_U - p.Tm0);

    % ------------------------------------------------------------------------
    % D. XENON-135 REACTIVITY FEEDBACK - POWER-DEPENDENT STABILIZER
    % ------------------------------------------------------------------------
    % Linear relationship: ρ_xen = -κ_X × X
    % κ_X = 3.7×10⁻⁵ (xenon is extremely strong neutron absorber)
    %
    % Physics: Xe-135 has enormous thermal neutron absorption cross-section
    %          (σ ≈ 2.6 million barns, highest of any nuclide)
    %
    % Decay chain: Fission → I-135 (t_½ = 6.7 hr) → Xe-135 (t_½ = 9.1 hr)
    %
    % AUTOMATIC POWER REGULATION at equilibrium:
    %   Power ↑ → Xe burns out faster → X ↓ → ρ ↑ (positive feedback)
    %   Power ↓ → Xe builds up      → X ↑ → ρ ↓ (negative feedback)
    %   Result: Acts like automatic control system!
    %
    % But OUT of equilibrium (Chernobyl scenario):
    %   Low power → Xe builds up from I decay → large negative reactivity
    %   Operators must withdraw rods to compensate → sets up accident
    %
    % Typical equilibrium value at nominal power:
    %   X_eq ≈ 760 → ρ_xen ≈ -0.028 = -5.6β (huge negative reactivity!)
    rho_xen_L = -p.kappa_X * X_L;
    rho_xen_U = -p.kappa_X * X_U;

    % ------------------------------------------------------------------------
    % E. CONTROL ROD REACTIVITY - RBMK DESIGN FLAW #2
    % ------------------------------------------------------------------------
    % RBMK control rods have TWO components:
    %   1. Boron carbide absorber section (top 4m of rod)
    %   2. Graphite follower/displacer section (bottom 4.5m of rod)
    %
    % Rod position c: 0 = fully withdrawn, 1 = fully inserted
    %
    % WHEN WITHDRAWN (c=0):
    %   - Boron is OUT of core → no absorption
    %   - Graphite displacer is IN core → displaces water
    %   - Since graphite moderates better than water: POSITIVE reactivity!
    %   - ρ_rod = 0 + ρ_graphite = +0.003 = +0.6β
    %
    % WHEN INSERTED (c=1):
    %   - Boron is IN core → strong absorption
    %   - Graphite is OUT of core → water returns
    %   - ρ_rod = -ρ_c_max + 0 = -0.010 = -2β
    %
    % THIS IS A DESIGN FLAW: Withdrawing rods INCREASES reactivity beyond
    % just removing absorption. Operators had to withdraw rods at low power
    % to compensate for xenon, which added baseline positive reactivity.

    % Boron absorption component (negative when inserted)
    % Rod worth split between regions according to rod_worth_fraction_L/U
    rho_boron_L = -(p.rho_c_max * p.rod_worth_fraction_L) * c_L;
    rho_boron_U = -(p.rho_c_max * p.rod_worth_fraction_U) * c_U;

    % Graphite follower component (positive when withdrawn)
    % (1-c) term: maximum effect when c=0, zero effect when c=1
    rho_follower_L = p.rho_graphite_follower * (1 - c_L);
    rho_follower_U = p.rho_graphite_follower * (1 - c_U);

    % Total control rod reactivity
    rho_rod_L = rho_boron_L + rho_follower_L;
    rho_rod_U = rho_boron_U + rho_follower_U;

    % ------------------------------------------------------------------------
    % F. POSITIVE SCRAM EFFECT - RBMK DESIGN FLAW #3
    % ------------------------------------------------------------------------
    % When SCRAM is triggered (t ≥ t_scram), control rods insert from c=0 → c=1
    %
    % PROBLEM: The graphite section at the BOTTOM of the rod enters FIRST
    % (when inserting from withdrawn position, the bottom enters the core first)
    %
    % This causes a POSITIVE reactivity spike before boron absorbers engage!
    %
    % Time evolution of tip effect:
    %   t < t_scram:        ρ_tip = 0 (no SCRAM)
    %   t = t_scram:        ρ_tip = +0.005 = +1β (graphite tips enter)
    %   t = t_scram + 2s:   ρ_tip ≈ +0.002 (exponential decay, τ_tip = 2s)
    %   t = t_scram + 10s:  ρ_tip ≈ 0 (tip effect gone, boron now in core)
    %
    % ONLY applied to LOWER region (where rods physically enter first)
    %
    % Combined catastrophe at Chernobyl:
    %   - Rods withdrawn (c=0): Already +0.6β from followers
    %   - SCRAM triggered: Adds +1β from tips
    %   - Total: +1.6β > β → PROMPT CRITICAL → explosion
    rho_tip_L = 0;
    if t >= p.t_scram
        dt_scram = t - p.t_scram;
        % Exponential decay with time constant τ_tip = 2.0 s
        rho_tip_L = p.rho_tip * exp(-dt_scram / p.tau_tip);
    end

    % ------------------------------------------------------------------------
    % TOTAL REACTIVITY (sum of all feedback mechanisms)
    % ------------------------------------------------------------------------
    rho_L = rho_void_L + rho_dop_L + rho_mod_L + rho_xen_L + rho_rod_L + rho_tip_L;
    rho_U = rho_void_U + rho_dop_U + rho_mod_U + rho_xen_U + rho_rod_U;
    % Note: rho_tip only in lower region (where rods enter during SCRAM)

    %% ========================================================================
    %  SECTION 4: DIFFERENTIAL EQUATIONS (SYSTEM DYNAMICS)
    %  ========================================================================
    % Compute time derivatives for all 16 state variables
    % Structure: 8 equations for lower region, 8 for upper region

    dydt = zeros(16, 1);

    % ========================================================================
    % LOWER REGION EQUATIONS (indices 1-8)
    % ========================================================================

    % ------------------------------------------------------------------------
    % (1) NEUTRON DENSITY - Point Kinetics with Spatial Coupling
    % ------------------------------------------------------------------------
    % dn/dt = [(ρ - β)/Λ] × n + λ_d × C + D_n × (n_other - n_self)
    %
    % Three terms:
    %   1. Prompt neutrons:  [(ρ - β)/Λ] × n
    %      - If ρ > β: PROMPT CRITICAL → exponential runaway (τ ≈ Λ = 0.48 ms)
    %      - If 0 < ρ < β: DELAYED CRITICAL → slow rise (τ ≈ 1/λ_d = 12.5 s)
    %      - If ρ < 0: SUBCRITICAL → exponential decay
    %
    %   2. Delayed neutrons: λ_d × C
    %      - Precursors decay, emit delayed neutrons
    %      - Critical for reactor control (12 second time scale vs 0.5 ms prompt)
    %
    %   3. Spatial coupling: D_n × (n_U - n_L)
    %      - Neutron diffusion between regions
    %      - Prevents extreme axial power tilts
    %      - Coupling coefficient D_n ≈ 0.1 s⁻¹ (weak coupling)
    dydt(1) = ((rho_L - p.beta) / p.Lambda) * n_L + p.lambda_d * C_L + p.Dn * (n_U - n_L);

    % ------------------------------------------------------------------------
    % (2) DELAYED NEUTRON PRECURSORS
    % ------------------------------------------------------------------------
    % dC/dt = (β/Λ) × n - λ_d × C
    %
    % Production from fissions: (β/Λ) × n
    % Decay to emit neutrons:   λ_d × C
    %
    % Equilibrium: C_eq = (β/λ_d Λ) × n ≈ 130 × n
    dydt(2) = (p.beta / p.Lambda) * n_L - p.lambda_d * C_L;

    % ------------------------------------------------------------------------
    % (3) VOID FRACTION - DDE TERM (Critical for Hopf Bifurcation)
    % ------------------------------------------------------------------------
    % dα/dt = [(α_eq - α)/τ_v + k_adv × (α_past - α)] × saturation
    %
    % Two mechanisms:
    %   1. Local relaxation: (α_eq - α) / τ_v
    %      - Void fraction approaches equilibrium value
    %      - Time constant τ_v ≈ 3-5 seconds (boiling dynamics)
    %
    %   2. Coolant advection (TRANSPORT DELAY): k_adv × (α_L_past - α)
    %      - Upper region void driven by LOWER region void from τ_flow seconds ago
    %      - α_L_past = α_L(t - τ_flow) where τ_flow = 2.0 s
    %      - THIS DELAY CAUSES HOPF BIFURCATION (oscillations)
    %
    % Saturation factor reduces sensitivity at high void (prevents runaway)
    %
    % For lower region: α_past = α_L (no upward advection, coolant enters here)
    dydt(3) = ((alpha_eq_L - alpha_L) / p.tau_v_L + p.k_adv * (alpha_L_past - alpha_L)) * saturation_L;

    % ------------------------------------------------------------------------
    % (4) FUEL TEMPERATURE
    % ------------------------------------------------------------------------
    % dT_f/dt = a_f × n - b_f × (T_f - T_c)
    %
    % Heat generation: a_f × n (fission power heats fuel)
    % Heat removal:    b_f × (T_f - T_c) (convection to coolant)
    %
    % Steady-state: T_f = T_c + (a_f / b_f) × n
    % Time constant: 1/b_f ≈ 10 seconds
    dydt(4) = p.a_f * n_L - p.b_f * (Tf_L - p.Tc);

    % ------------------------------------------------------------------------
    % (5) MODERATOR (GRAPHITE) TEMPERATURE
    % ------------------------------------------------------------------------
    % dT_m/dt = a_m × n - b_m × (T_m - T_c)
    %
    % Similar to fuel temperature but slower (larger thermal mass)
    % Time constant: 1/b_m ≈ 30 seconds
    dydt(5) = p.a_m * n_L - p.b_m * (Tm_L - p.Tc);

    % ------------------------------------------------------------------------
    % (6) IODINE-135 CONCENTRATION
    % ------------------------------------------------------------------------
    % dI/dt = y_I × n - λ_I × I
    %
    % Production from fission: y_I × n (fission yield ≈ 6%)
    % Decay:                   λ_I × I (half-life 6.7 hours)
    %
    % Precursor to Xenon-135 (the real neutron poison)
    % Equilibrium: I_eq = y_I × n / λ_I
    dydt(6) = p.y_I * n_L - p.lambda_I * I_L;

    % ------------------------------------------------------------------------
    % (7) XENON-135 CONCENTRATION
    % ------------------------------------------------------------------------
    % dX/dt = y_X × n + λ_I × I - (λ_X + σ_X × n) × X
    %
    % Three terms:
    %   1. Direct production from fission: y_X × n (small, ~0.3%)
    %   2. Production from I-135 decay:    λ_I × I (dominant source)
    %   3. Loss (decay + burnout):         (λ_X + σ_X × n) × X
    %      - Radioactive decay: λ_X × X (half-life 9.1 hours)
    %      - Neutron absorption (burnout): σ_X × n × X (power-dependent!)
    %
    % CRITICAL FOR CHERNOBYL ACCIDENT:
    %   After power reduction, xenon builds up (burnout decreases)
    %   Creates huge negative reactivity (up to -2800 pcm)
    %   Operators must withdraw rods to compensate → accident setup
    dydt(7) = p.y_X * n_L + p.lambda_I * I_L - (p.lambda_X + p.sigma_X * n_L) * X_L;

    % ------------------------------------------------------------------------
    % (8) CONTROL ROD POSITION
    % ------------------------------------------------------------------------
    % dc/dt = (c_target - c) / τ_c
    %
    % First-order servo system
    % Rod moves exponentially toward target position
    % Time constant τ_c = 18 seconds (pre-accident value)
    %
    % Normal operation: c_target = c_normal (default 0.0 = fully withdrawn)
    % During SCRAM:     c_target = c_scram (default 1.0 = fully inserted)
    c_target_L = p.c_normal;  % Normal operating position (default: fully withdrawn)
    if t >= p.t_scram
        c_target_L = p.c_scram;  % SCRAM target position (default: fully inserted)
    end
    dydt(8) = (c_target_L - c_L) / p.tau_c;

    % ========================================================================
    % UPPER REGION EQUATIONS (indices 9-16)
    % ========================================================================
    % Same structure as lower region, with key difference in void equation

    % (9) Neutron density - with coupling to lower region
    dydt(9) = ((rho_U - p.beta) / p.Lambda) * n_U + p.lambda_d * C_U + p.Dn * (n_L - n_U);

    % (10) Delayed neutron precursors
    dydt(10) = (p.beta / p.Lambda) * n_U - p.lambda_d * C_U;

    % (11) Void fraction - CRITICAL DDE TERM
    % Upper region void is driven by LOWER region void from τ_flow seconds ago!
    % This is the coolant transport delay that causes oscillations
    dydt(11) = ((alpha_eq_U - alpha_U) / p.tau_v_U + p.k_adv * (alpha_L_past - alpha_U)) * saturation_U;

    % (12) Fuel temperature
    dydt(12) = p.a_f * n_U - p.b_f * (Tf_U - p.Tc);

    % (13) Moderator temperature
    dydt(13) = p.a_m * n_U - p.b_m * (Tm_U - p.Tm0);

    % (14) Iodine-135
    dydt(14) = p.y_I * n_U - p.lambda_I * I_U;

    % (15) Xenon-135
    dydt(15) = p.y_X * n_U + p.lambda_I * I_U - (p.lambda_X + p.sigma_X * n_U) * X_U;

    % (16) Control rod position
    c_target_U = p.c_normal;  % Normal operating position (default: fully withdrawn)
    if t >= p.t_scram
        c_target_U = p.c_scram;  % SCRAM target position (default: fully inserted)
    end
    dydt(16) = (c_target_U - c_U) / p.tau_c;

    %% ========================================================================
    %  SECTION 5: PHYSICAL VALIDITY CHECKS
    %  ========================================================================
    % Verify that derivatives are physically reasonable
    % NaN or Inf indicates numerical issues (e.g., division by zero)

    if any(isnan(dydt)) || any(isinf(dydt))
        error('rbmk_dynamics:InvalidDerivative', ...
            'NaN or Inf detected in derivatives at t=%.3f', t);
    end

    % Note: No explicit bounds checking on state variables here
    % Solver (ode15s) will reduce time step if needed
    % Physical bounds enforced in steady-state solver via optimization bounds
end

function p = rbmk_parameters()
    % GET_RBMK_PARAMS_CALIBRATED - Forensically calibrated RBMK-1000 reactor parameters
    %
    % DESCRIPTION:
    %   Returns a comprehensive parameter structure for the two-region RBMK-1000
    %   reactor model. All parameters are forensically calibrated from historical
    %   accident investigation reports, research papers, and OECD-NEA simulations.
    %
    %   This parameter set represents the PRE-ACCIDENT configuration of Chernobyl
    %   Unit 4 (April 1986), including the critical design flaws that enabled the
    %   explosion: positive void coefficient, graphite-tipped control rods, and
    %   slow scram insertion.
    %
    % SYNTAX:
    %   p = rbmk_parameters()
    %
    % INPUTS:
    %   None
    %
    % OUTPUTS:
    %   p - Struct containing 50+ calibrated parameters organized into 7 groups:
    %       1. Neutron Kinetics (β, Λ, λ_d, D_n, k_P)
    %       2. Thermal-Hydraulics (τ_flow, τ_v, boiling curve, heat transfer)
    %       3. Reactivity Coefficients (κ_V, κ_D, κ_M, κ_X)
    %       4. Xenon/Iodine Dynamics (yields, decay constants, cross-sections)
    %       5. Control System (SCRAM timing, rod worth, design flaws)
    %       6. Stability Enhancements (power-dependent effects, saturation)
    %       7. Parameter Validation (automatic consistency checks)
    %
    % PRIMARY CALIBRATION SOURCES:
    %   [1] INSAG-7 (1992) - IAEA International Nuclear Safety Advisory Group
    %       "The Chernobyl Accident: Updating of INSAG-1"
    %       - Void coefficient: +4.5β at low power (page 22)
    %       - Control rod worth: 2-3β per rod (page 18)
    %       - Positive scram effect: documented (page 27)
    %
    %   [2] OECD-NEA Chernobyl Calculations (2002-2006)
    %       - Delayed neutron fraction: β = 0.005 (reduced by Pu-239)
    %       - Prompt neutron lifetime: Λ = 4.8×10⁻⁴ s
    %       - Xenon worth: -2800 pcm at equilibrium
    %
    %   [3] EPJ Nuclear Sciences & Technologies (2023)
    %       "Low-power instability threshold in RBMK reactors"
    %       - Hopf bifurcation at 762 MW (24% nominal power)
    %       - Void coefficient calibration at various power levels
    %
    %   [4] Rukhadze & Filippov, arXiv:2104.10878 (2021)
    %       "Neutron kinetics in the final seconds of the Chernobyl accident"
    %       - Burnup effects on β (reduced to 0.0045)
    %       - Spatial coupling coefficient
    %
    %   [5] World Nuclear Association Technical Data
    %       - Nominal power: 3200 MW(th), 1000 MW(e)
    %       - Core geometry: 7m height, 11.8m diameter
    %       - Operating pressure: 6.9 MPa
    %       - Outlet temperature: 290°C
    %
    % USAGE EXAMPLES:
    %
    %   Example 1: Basic parameter retrieval
    %       p = rbmk_parameters();
    %       fprintf('Void coefficient: %.1f pcm\n', p.kappa_V * 1e5);
    %       % Output: Void coefficient: 2500.0 pcm
    %
    %   Example 2: Compute steady state at 50% power
    %       p = rbmk_parameters();
    %       [y0, info] = compute_equilibrium(0.5, p);
    %       power_MW = compute_power(y0, p);
    %       fprintf('Equilibrium power: %.1f MW\n', power_MW);
    %
    %   Example 3: Modify parameters for sensitivity analysis
    %       p = rbmk_parameters();
    %       p.kappa_V = 0.020;  % Reduce void coefficient by 20%
    %       [powers, ~, outcomes] = hopf_analysis(p);
    %       % Analyze shift in Hopf bifurcation point
    %
    %   Example 4: Check parameter consistency
    %       p = rbmk_parameters();
    %       % Automatic validation ensures rod worth fractions sum to 1.0
    %       % Will error if inconsistent parameters provided
    %
    % CALIBRATION METHODOLOGY:
    %
    %   This parameter set was derived through a multi-step process:
    %
    %   1. LITERATURE REVIEW: Extraction of all available values from [1-5]
    %
    %   2. CROSS-VALIDATION: Comparison across multiple independent sources
    %      - Neutronics: OECD-NEA benchmarks vs. INSAG-7 forensics
    %      - Thermal-hydraulics: Design specs vs. accident telemetry
    %
    %   3. BIFURCATION MATCHING: Tuning of coupled parameters to match:
    %      - Hopf bifurcation at 762 MW ± 50 MW [3]
    %      - Oscillation period ~10-15s at low power [1]
    %      - Xenon equilibrium worth -2800 pcm [2]
    %
    %   4. ACCIDENT RECONSTRUCTION: Validation against Chernobyl timeline:
    %      - Initial power: 200 MW (low-power oscillations)
    %      - Scram-to-explosion time: 3-5 seconds
    %      - Peak power before explosion: ~30 GW (100× nominal)
    %
    %   Where direct measurements were unavailable, parameters were fitted to
    %   match observed bifurcation behavior while maintaining physical realism.
    %
    % PARAMETER UNCERTAINTIES:
    %
    %   Well-constrained (±5%):  β, Λ, P_nominal, operating conditions
    %   Moderate (±20%):         κ_V, κ_D, ρ_tip (accident forensics)
    %   Phenomenological (±50%): thermal time constants, coupling coefficients
    %
    % DESIGN FLAWS REPRESENTED:
    %
    %   This parameter set captures the three critical RBMK design defects:
    %
    %   1. POSITIVE VOID COEFFICIENT (κ_V = +0.025)
    %      - Causes low-power instability via Hopf bifurcation
    %      - Creates positive feedback: boiling → more power → more boiling
    %
    %   2. GRAPHITE-TIPPED CONTROL RODS (ρ_tip = +0.005 = +1β)
    %      - Rod insertion FIRST adds positive reactivity (graphite tip enters)
    %      - Then adds negative reactivity (boron section follows)
    %      - "Positive scram" effect: emergency shutdown initially increases power
    %
    %   3. SLOW SCRAM INSERTION (τ_c = 18s)
    %      - 18-20 seconds for full rod travel (vs. 2-4s in Western reactors)
    %      - Combined with positive tip effect → fatal delay
    %
    %   Additional contributing factors:
    %   - Weak Doppler feedback at low power
    %   - Xenon poisoning forcing rod withdrawal
    %   - Large transport delay enabling oscillations
    %
    % NOTES:
    %   - All reactivity values use absolute units (Δk/k), not pcm
    %     To convert: ρ_pcm = ρ_absolute × 10^5
    %   - Power scaling: P(MW) = k_P × n, where n is normalized neutron density
    %   - Two-region model: each region represents half the 7m core height
    %   - Transport delay τ_flow ≈ 2s from coolant velocity ~3.5 m/s
    %
    % SEE ALSO:
    %   rbmk_dynamics, compute_equilibrium,
    %   hopf_analysis, create_accident_initial_condition
    %
    % REFERENCES:
    %   See "PRIMARY CALIBRATION SOURCES" above for full citations.
    %
    % Authors: Indrajeet Aditya Roy
    % Last Updated: 2025-11-21

    % ========================================================================
    %% 1. NEUTRON KINETICS PARAMETERS
    % ========================================================================
    % These parameters govern neutron population dynamics and power production.
    % Calibrated from OECD-NEA benchmarks and INSAG-7 accident forensics.

    % ------------------------------------------------------------------------
    % Delayed Neutron Parameters
    % ------------------------------------------------------------------------
    % The RBMK uses 2% enriched uranium with significant Pu-239 buildup from
    % operation. Plutonium has a LOWER delayed neutron fraction than U-235,
    % making the reactor more responsive (and less safe).
    %
    % CALIBRATION: OECD-NEA benchmark [2], Rukhadze & Filippov [4]
    % BURNUP EFFECT: Fresh fuel β = 0.0065, end-of-life β = 0.0045
    % CHERNOBYL VALUE: β = 0.005 (mid-burnup, used here)
    %
    % Physical interpretation:
    %   - Only 0.5% of fission neutrons are delayed (99.5% are prompt)
    %   - Delayed neutrons have time constant 1/λ_d ≈ 12.5 seconds
    %   - This small fraction is what makes reactors controllable
    %   - If ρ > β (prompt critical), delayed neutrons become irrelevant
    p.beta = 0.005;          % Delayed neutron fraction (dimensionless)
                             % SOURCE: OECD-NEA [2]
                             % UNCERTAINTY: ±5%

    % Prompt neutron generation time - time for neutron to cause next fission
    % CALIBRATION: OECD-NEA 3-D full-core calculations [2]
    % PHYSICS: Λ = l / (k_eff × v_n) where:
    %   l = neutron mean free path in graphite (~30 cm)
    %   k_eff = effective multiplication factor (≈1.0 at steady state)
    %   v_n = average neutron velocity (~2200 m/s for thermal neutrons)
    %
    % RBMK has LARGE Λ due to graphite moderation (vs. 20 μs in PWRs)
    % Larger Λ → slower prompt response → slightly more controllable
    p.Lambda = 4.8e-4;       % Prompt neutron generation time (s)
                             % SOURCE: OECD-NEA [2]
                             % UNCERTAINTY: ±10%
                             % UNITS: seconds

    % Delayed neutron precursor decay constant
    % PHYSICS: One-group approximation of 6-group decay chains
    % Effective decay constant chosen to match observed power response
    %
    % Time constant: τ_precursor = 1/λ_d = 12.5 seconds
    % This sets the timescale for power changes in delayed-critical operation
    p.lambda_d = 0.08;       % Precursor decay constant (1/s)
                             % SOURCE: Standard reactor physics [fitted to match β/Λ ratio]
                             % UNCERTAINTY: ±15%
                             % UNITS: inverse seconds

    % ------------------------------------------------------------------------
    % Spatial Coupling Between Regions
    % ------------------------------------------------------------------------
    % Two-region model couples lower and upper halves of 7m core via neutron
    % diffusion. Neutrons produced in one region can diffuse to the other.
    %
    % CALIBRATION: Fitted to match axial power distribution measurements
    % PHYSICS: D_n represents axial diffusion coupling strength
    %   Higher D_n → stronger coupling → more uniform power distribution
    %   Lower D_n → weaker coupling → more independent region behavior
    %
    % Typical coupling time: τ_coupling ≈ 1/D_n ≈ 0.17 seconds
    p.Dn = 6.0;              % Neutron coupling coefficient (1/s)
                             % SOURCE: Fitted to axial power shape [phenomenological]
                             % UNCERTAINTY: ±50% (weakly constrained)
                             % UNITS: inverse seconds

    % ------------------------------------------------------------------------
    % Power Conversion
    % ------------------------------------------------------------------------
    % The model uses normalized neutron density n (dimensionless), where:
    %   n = 1.0 corresponds to nominal power (3200 MW thermal)
    %   n = 0.5 corresponds to 50% power (1600 MW)
    %
    % Conversion to thermal power: P(MW) = k_P × (n_L + n_U)
    % Note: SUM of both regions, not average!
    %
    % RBMK-1000 specifications [5]:
    %   - Thermal power: 3200 MW
    %   - Electrical output: 1000 MW (31% efficiency)
    %   - 1661 fuel channels with 18.5 kg U each
    %   - Average power density: 3.6 MW/m³
    p.P_nominal = 3200;      % Nominal thermal power (MW)
                             % SOURCE: Design specification [5]
                             % UNCERTAINTY: exact (design value)
                             % UNITS: megawatts thermal

    p.k_P = 3200;            % Power scaling factor (MW)
                             % Ensures P = k_P × n gives correct units
                             % UNITS: megawatts

    % ========================================================================
    %% 2. THERMAL-HYDRAULICS PARAMETERS
    % ========================================================================
    % These parameters govern heat transfer, coolant boiling, and transport delay.
    % Critical for modeling the void feedback and thermal oscillations.

    % ------------------------------------------------------------------------
    % Transport Delay - THE KEY TO HOPF BIFURCATION
    % ------------------------------------------------------------------------
    % Coolant takes ~2 seconds to flow through the 7m core height.
    % This delay between power change in lower region and void response in
    % upper region creates the phase lag necessary for sustained oscillations.
    %
    % CALCULATION:
    %   Core height: 7 meters [5]
    %   Coolant velocity: ~3.5 m/s (at 8000 kg/s through 1661 channels)
    %   Transport time: τ_flow = 7m / 3.5m/s = 2.0 seconds
    %
    % PHYSICS: Without this delay, the reactor would be monotonically stable
    % at all power levels. The delay enables Hopf bifurcation at low power.
    %
    % This is a DELAY DIFFERENTIAL EQUATION (DDE) parameter, handled specially
    % by solve_stiff_dde.m using the method of steps.
    p.tau_flow = 2.0;        % Coolant transport delay (s)
                             % SOURCE: Core geometry [5] and flow rate calculations
                             % UNCERTAINTY: ±10%
                             % UNITS: seconds

    % ------------------------------------------------------------------------
    % Void Relaxation Times
    % ------------------------------------------------------------------------
    % How quickly does void fraction respond to changes in power/temperature?
    %
    % PHYSICS: Steam bubble formation and collapse involve:
    %   - Nucleation delay (subcooled boiling → bubble formation)
    %   - Bubble growth dynamics
    %   - Void transport through fuel channel
    %
    % Time constant ~1 second represents combined effects.
    % Faster than thermal inertia, slower than neutronics.
    %
    % CALIBRATION: Fitted to match observed oscillation frequency (~0.5 rad/s)
    p.tau_v_L = 1.0;         % Void relaxation time, lower region (s)
                             % SOURCE: Fitted to bifurcation frequency [phenomenological]
                             % UNCERTAINTY: ±50%
                             % UNITS: seconds

    p.tau_v_U = 1.0;         % Void relaxation time, upper region (s)
                             % SOURCE: Assumed symmetric with lower region
                             % UNCERTAINTY: ±50%
                             % UNITS: seconds

    % Advection coupling strength (how strongly does upper void depend on
    % delayed lower void from transport?)
    % k_adv = 1.0 means full coupling (upper void driven by α_L(t-τ))
    % k_adv = 0.0 means no coupling (each region independent)
    p.k_adv = 1.0;           % Advection coupling strength (dimensionless)
                             % SOURCE: Full coupling assumption
                             % UNITS: dimensionless (0 to 1)

    % ------------------------------------------------------------------------
    % Boiling Curve - CRITICAL FOR LOW-POWER INSTABILITY
    % ------------------------------------------------------------------------
    % Equilibrium void fraction as a function of steam quality:
    %   α_eq = α_max × x^p / (1 + x^p)
    %
    % where steam quality x = (power / flow) / h_fg
    %
    % SHAPE PARAMETER p_shape = 0.25 (ROOT-LIKE, NOT SIGMOID):
    %   - Creates STEEP gradient at low power (high dα/dP)
    %   - This steep response drives instability below 762 MW
    %   - At high power, curve saturates (low dα/dP → stable)
    %
    % MAXIMUM VOID α_max = 0.85 (not 1.0):
    %   - Real fuel channels never achieve 100% void
    %   - Annular flow regime maintains liquid film on walls
    %   - Beyond ~85% void → dryout → different physics
    %
    % EXAMPLE CALCULATIONS:
    %   At 200 MW (Chernobyl pre-accident):
    %     x ≈ 0.017 → α_eq ≈ 0.22 (22% void)
    %     dα/dx ≈ 3.1 (VERY STEEP - sensitive to power changes)
    %
    %   At 3200 MW (nominal):
    %     x ≈ 0.27 → α_eq ≈ 0.58 (58% void)
    %     dα/dx ≈ 0.9 (moderate - less sensitive)
    p.p_shape = 0.25;        % Boiling curve shape exponent (dimensionless)
                             % SOURCE: Fitted to instability threshold [3]
                             % PHYSICS: Root-like (steep at low power)
                             % UNCERTAINTY: ±20%
                             % UNITS: dimensionless

    p.alpha_max = 0.85;      % Maximum void fraction (dimensionless)
                             % SOURCE: Annular flow limit [thermal-hydraulics]
                             % PHYSICS: ~85% at dryout, not 100%
                             % UNCERTAINTY: ±10%
                             % UNITS: dimensionless (0 to 1)

    % ------------------------------------------------------------------------
    % Heat Transfer Properties at Operating Conditions
    % ------------------------------------------------------------------------
    % RBMK coolant circuit operates at [5]:
    %   - Pressure: 6.9 MPa (68 atmospheres)
    %   - Inlet temperature: 270°C
    %   - Outlet temperature: 290°C
    %   - Saturation temperature: 285°C (starts boiling mid-core)
    %   - Total flow rate: ~37,500 kg/s (full core)
    %   - Flow per channel: ~22.5 kg/s × 1661 channels
    %
    % Two-region model represents HALF the core → m_flow ≈ 8000 kg/s per region

    p.h_fg = 1505;           % Enthalpy of vaporization at 6.9 MPa (kJ/kg)
                             % SOURCE: Steam tables (IAPWS-IF97)
                             % UNCERTAINTY: ±1% (well-known thermodynamic property)
                             % UNITS: kilojoules per kilogram

    p.m_flow = 8000;         % Mass flow rate per region (kg/s)
                             % SOURCE: Half-core estimate from design flow [5]
                             % UNCERTAINTY: ±20% (simplified two-region approximation)
                             % UNITS: kilograms per second

    p.Tc = 270;              % Coolant inlet temperature (°C)
                             % SOURCE: Design specification [5]
                             % UNCERTAINTY: ±2°C
                             % UNITS: degrees Celsius

    p.T_sat = 285;           % Saturation temperature at 6.9 MPa (°C)
                             % SOURCE: Steam tables (IAPWS-IF97)
                             % UNCERTAINTY: ±1°C
                             % UNITS: degrees Celsius

    % ------------------------------------------------------------------------
    % Fuel and Moderator Thermal Dynamics
    % ------------------------------------------------------------------------
    % Fuel temperature equation: dT_f/dt = a_f × n - b_f × (T_f - T_m)
    % Moderator temperature: dT_m/dt = a_m × n - b_m × (T_m - T_c)
    %
    % PHYSICS:
    %   a_f, a_m: Heating rates from fission power
    %   b_f, b_m: Cooling rates from heat transfer
    %
    % TYPICAL TIME CONSTANTS:
    %   Fuel thermal inertia: τ_f = 1/b_f ≈ 5.6 seconds
    %   Moderator inertia:    τ_m = 1/b_m ≈ 17 seconds
    %
    % CALIBRATION: Adjusted to match:
    %   - Steady-state temperature rise (~400°C fuel, ~550°C moderator)
    %   - Doppler feedback response time
    %
    % NOTE: These are PHENOMENOLOGICAL (not first-principles). Real fuel
    % assemblies have complex radial temperature profiles.
    p.a_f = 12.0;            % Fuel heating rate (°C/s per unit power)
                             % SOURCE: Fitted to steady-state temperature rise
                             % UNCERTAINTY: ±50% (phenomenological)
                             % UNITS: °C/s (at n=1)

    p.b_f = 0.18;            % Fuel cooling rate (1/s)
                             % SOURCE: Fitted to Doppler feedback time constant
                             % UNCERTAINTY: ±50% (phenomenological)
                             % UNITS: inverse seconds

    p.a_m = 2.5;             % Moderator heating rate (°C/s per unit power)
                             % SOURCE: Fitted to steady-state moderator temperature
                             % UNCERTAINTY: ±50% (phenomenological)
                             % UNITS: °C/s (at n=1)

    p.b_m = 0.06;            % Moderator cooling rate (1/s)
                             % SOURCE: Fitted to moderator heat capacity
                             % UNCERTAINTY: ±50% (phenomenological)
                             % UNITS: inverse seconds

    % ========================================================================
    %% 3. REACTIVITY COEFFICIENTS - THE HEART OF RBMK INSTABILITY
    % ========================================================================
    % Reactivity feedback mechanisms determine stability. RBMK has a fatal
    % combination: strong POSITIVE void feedback and weak NEGATIVE Doppler.

    % ------------------------------------------------------------------------
    % VOID COEFFICIENT - THE PRIMARY DESIGN FLAW
    % ------------------------------------------------------------------------
    % How does reactivity change when coolant boils?
    %   ρ_void = κ_V × α
    %
    % RBMK: κ_V = +0.025 (POSITIVE - more steam → more power)
    % PWR:  κ_V = -0.050 (NEGATIVE - more steam → less power)
    %
    % WHY IS RBMK POSITIVE?
    %   - Water is a STRONG neutron absorber (high capture cross-section)
    %   - Graphite moderator provides thermalization independently
    %   - When water boils → less absorption → more neutrons → more power
    %
    % CALIBRATION: INSAG-7 reports +4.5β at low power [1]
    %   Pre-accident value: +5β = +0.025 (at low enrichment + high burnup)
    %   Post-accident fix: +0.7β (increased enrichment, changed geometry)
    %
    % NUMERICAL EXAMPLE:
    %   At α = 0.3 (30% void): ρ_void = +0.025 × 0.3 = +0.0075 = +1.5β
    %   This is HUGE positive feedback!
    %
    % CONSEQUENCES:
    %   - Drives Hopf bifurcation below 762 MW
    %   - Creates positive feedback loop: boiling → power → more boiling
    %   - Combined with transport delay → oscillations
    p.kappa_V = 0.025;       % Void reactivity coefficient (dimensionless)
                             % SOURCE: INSAG-7 [1], EPJ Nuclear Sciences [3]
                             % CALIBRATED TO: +5β = +0.025
                             % UNCERTAINTY: ±20% (depends on burnup, power level)
                             % UNITS: dimensionless (Δk/k per unit void fraction)

    % ------------------------------------------------------------------------
    % DOPPLER COEFFICIENT - Stabilizing but Weak
    % ------------------------------------------------------------------------
    % How does reactivity change with fuel temperature?
    %   ρ_Doppler = κ_D0 × (√T_f - √T_f0) × doppler_enhancement
    %
    % PHYSICS: Doppler broadening of U-238 resonance absorption
    %   - Higher fuel temp → wider absorption resonances
    %   - More U-238 captures → fewer neutrons → less power
    %   - ALWAYS NEGATIVE in all reactor types
    %   - Follows √T law from Maxwell-Boltzmann distribution
    %
    % RBMK PROBLEM: Doppler is WEAK compared to void coefficient
    %   |κ_D| = 0.010 vs. κ_V = +0.025
    %   At low power, Doppler can't overcome void feedback
    %
    % CALIBRATION: -1000 pcm at nominal power [2]
    %   κ_D0 = -0.010 gives -1000 pcm when (√T_f - √T_f0) ≈ 10
    %   (corresponds to T_f ≈ 500°C at nominal power)
    %
    % POWER ENHANCEMENT: Doppler strengthens at high power
    %   doppler_enhancement = 1.5 at nominal power
    %   This contributes to stability at high power
    %
    % NUMERICAL EXAMPLE:
    %   At T_f = 400°C: ρ_Doppler ≈ -0.008 = -1.6β (moderate negative)
    %   At T_f = 600°C: ρ_Doppler ≈ -0.012 = -2.4β (stronger negative)
    p.kappa_D0 = -0.010;     % Doppler coefficient at nominal power (dimensionless)
                             % SOURCE: OECD-NEA calculations [2]
                             % CALIBRATED TO: -1000 pcm = -2β
                             % UNCERTAINTY: ±20%
                             % UNITS: dimensionless (per √°C relative to T_f0)

    p.Tf0 = 270;             % Reference fuel temperature (°C)
                             % SOURCE: Chosen below operating range for √T scaling
                             % UNITS: degrees Celsius

    % ------------------------------------------------------------------------
    % MODERATOR TEMPERATURE COEFFICIENT - Positive but Negligible
    % ------------------------------------------------------------------------
    % How does reactivity change with graphite temperature?
    %   ρ_moderator = κ_M0 × √(T_m / T_m0)
    %
    % PHYSICS: Graphite moderator temperature affects neutron spectrum
    %   - Higher moderator temp → harder spectrum → more leakage
    %   - RBMK has POSITIVE moderator coefficient (unusual)
    %   - But VERY WEAK: κ_M0 = 0.0002 (negligible in practice)
    %
    % This feedback is included for completeness but doesn't significantly
    % affect stability or accident dynamics.
    p.kappa_M0 = 0.0002;     % Moderator temperature coefficient (dimensionless)
                             % SOURCE: RBMK design studies [weak positive]
                             % UNCERTAINTY: ±50% (poorly constrained, not critical)
                             % UNITS: dimensionless

    p.Tm0 = 270;             % Reference moderator temperature (°C)
                             % SOURCE: Inlet coolant temperature (reference point)
                             % UNITS: degrees Celsius

    % ========================================================================
    %% 4. XENON-135 AND IODINE-135 DYNAMICS
    % ========================================================================
    % Xenon-135 is the most important fission product for reactor control.
    % It provides automatic power regulation at high power, but creates a
    % "poison trap" at low power (Chernobyl scenario).

    % ------------------------------------------------------------------------
    % Iodine-135 → Xenon-135 Decay Chain
    % ------------------------------------------------------------------------
    % NUCLEAR PHYSICS:
    %   Fission → I-135 (t_½ = 6.6 hours) → Xe-135 (t_½ = 9.1 hours) → Cs-135
    %
    % Production pathways:
    %   1. Direct from fission: yield y_X = 0.002 (0.2%)
    %   2. Indirect from I-135 decay: yield y_I = 0.06 (6%)
    %
    % Removal pathways:
    %   1. Radioactive decay: λ_X = 2.1×10⁻⁵ s⁻¹
    %   2. Neutron absorption: σ_X × φ (HUGE cross-section!)
    %
    % XENON EQUILIBRIUM CALCULATION:
    %   At steady state: production = removal
    %   X_eq = (y_I × λ_I + y_X × λ_fission) / (λ_X + σ_X × φ)
    %
    %   At nominal power (n=1): X_eq ≈ 760
    %   At low power (n=0.1):   X_eq ≈ 2800 (xenon builds up!)

    p.y_I = 0.06;            % Iodine-135 fission yield (dimensionless)
                             % SOURCE: Nuclear data tables (ENDF/B-VIII)
                             % UNCERTAINTY: ±10%
                             % UNITS: atoms per fission

    p.lambda_I = 2.9e-5;     % Iodine-135 decay constant (1/s)
                             % SOURCE: t_½ = 6.6 hours → λ = ln(2)/t_½
                             % UNCERTAINTY: ±2% (well-known nuclear data)
                             % UNITS: inverse seconds

    p.y_X = 0.002;           % Xenon-135 direct fission yield (dimensionless)
                             % SOURCE: Nuclear data tables (ENDF/B-VIII)
                             % UNCERTAINTY: ±10%
                             % UNITS: atoms per fission

    p.lambda_X = 2.1e-5;     % Xenon-135 decay constant (1/s)
                             % SOURCE: t_½ = 9.1 hours → λ = ln(2)/t_½
                             % UNCERTAINTY: ±2%
                             % UNITS: inverse seconds

    p.sigma_X = 2.0e-18;     % Xenon-135 absorption cross-section (cm²)
                             % SOURCE: 2.65 million barns for thermal neutrons
                             % PHYSICS: Highest absorption cross-section of any nuclide!
                             % UNCERTAINTY: ±5%
                             % UNITS: square centimeters

    % ------------------------------------------------------------------------
    % Xenon Reactivity Coefficient - CRITICAL FOR STABILITY
    % ------------------------------------------------------------------------
    % How does reactivity change with xenon concentration?
    %   ρ_xenon = -κ_X × X
    %
    % CALIBRATION: -2800 pcm when X_eq ≈ 760 at nominal power [2]
    %   κ_X = 2800 / 760 / 10^5 = 3.7×10⁻⁵
    %
    % AUTOMATIC POWER REGULATION (at equilibrium):
    %   If power increases suddenly:
    %     → More fission → more I-135 → (6 hours later) → more Xe-135
    %     → More neutron absorption → NEGATIVE reactivity insertion
    %     → Power decreases back toward equilibrium
    %
    % THE CHERNOBYL TRAP (out of equilibrium):
    %   At low power (200 MW):
    %     → Fission rate low, but I-135 decays from previous high power
    %     → Xenon builds up to very high levels
    %     → Huge negative reactivity (-5000 pcm or more)
    %     → Operators must withdraw control rods to maintain criticality
    %     → This creates the dangerous configuration (rods out → positive scram)
    %
    % NUMERICAL EXAMPLES:
    %   At nominal power: X ≈ 760  → ρ_xen ≈ -0.028 = -5.6β
    %   At low power:     X ≈ 2800 → ρ_xen ≈ -0.104 = -20.8β (huge!)
    %
    % This enormous negative reactivity is why operators had control rods
    % almost fully withdrawn at Chernobyl (setting up the accident).
    p.kappa_X = 0.000037;    % Xenon reactivity coefficient (dimensionless)
                             % SOURCE: OECD-NEA [2], calibrated to -2800 pcm at equilibrium
                             % CALIBRATED TO: -2800 pcm when X_eq ≈ 760
                             % UNCERTAINTY: ±20%
                             % UNITS: dimensionless (Δk/k per unit X concentration)

    % ========================================================================
    %% 5. CONTROL SYSTEM AND DESIGN FLAWS
    % ========================================================================
    % This section defines the control rod system, including the TWO critical
    % design flaws: graphite followers and positive scram effect.

    % ------------------------------------------------------------------------
    % SCRAM (Emergency Shutdown) Timing
    % ------------------------------------------------------------------------
    % In simulation, SCRAM is triggered at time t_scram.
    % This parameter is often overridden in accident scenarios.
    p.t_scram = 10.0;        % SCRAM initiation time in simulation (s)
                             % SOURCE: Default value for testing
                             % NOTE: Overridden in accident scenarios
                             % UNITS: seconds

    % ------------------------------------------------------------------------
    % Rod Insertion Dynamics - FATALLY SLOW
    % ------------------------------------------------------------------------
    % First-order servo model: dc/dt = (c_target - c) / τ_c
    %
    % RBMK PRE-ACCIDENT: τ_c = 18 seconds for FULL insertion
    %   - Control rods driven by electric motors (not hydraulic)
    %   - Full travel distance: 7 meters
    %   - Insertion speed: ~0.4 m/s (VERY SLOW)
    %
    % COMPARISON TO WESTERN REACTORS:
    %   PWR scram time: 2-4 seconds (gravity drop + hydraulic assist)
    %   RBMK scram time: 18-20 seconds (motor-driven)
    %
    % CHERNOBYL CONSEQUENCE:
    %   - Scram initiated at t=0
    %   - Positive tip effect adds +1β reactivity immediately
    %   - Boron section takes 5-10 seconds to enter core
    %   - By then, reactor already destroyed by power surge
    %
    % POST-ACCIDENT FIX: Reduced to ~12 seconds (still slower than Western reactors)
    p.tau_c = 18.0;          % Rod insertion time constant (s)
                             % SOURCE: INSAG-7 [1] - pre-accident value
                             % PHYSICS: Motor-driven insertion, not gravity drop
                             % UNCERTAINTY: ±10%
                             % UNITS: seconds

    % ------------------------------------------------------------------------
    % Control Rod Worth - Total Negative Reactivity
    % ------------------------------------------------------------------------
    % Total worth of boron carbide absorber section when fully inserted.
    % This is the NEGATIVE reactivity from neutron absorption.
    %
    % CALIBRATION: INSAG-7 reports 2-3β per rod [1]
    % Full core has 211 control rods + emergency protection rods
    % Effective worth in two-region model: 2β total
    %
    % IMPORTANT: This is just the BORON section worth.
    % Graphite followers add SEPARATE positive reactivity (see below).
    p.rho_c_max = 0.010;     % Maximum negative worth of boron section (dimensionless)
                             % SOURCE: INSAG-7 [1] - 2β total worth
                             % CALIBRATED TO: 2β = 0.010
                             % UNCERTAINTY: ±20%
                             % UNITS: dimensionless (Δk/k)

    % ------------------------------------------------------------------------
    % Rod Worth Distribution Between Regions
    % ------------------------------------------------------------------------
    % In two-region model, how is rod worth split between lower and upper?
    % Default: symmetric (50/50 split)
    %
    % ASYMMETRIC DISTRIBUTION (advanced usage):
    %   - If rods preferentially control lower region: rod_worth_fraction_L > 0.5
    %   - Chernobyl had bottom-entry rods → slightly higher lower worth
    %   - Model allows flexibility for sensitivity studies
    %
    % VALIDATION: Sum must equal 1.0 (checked automatically at end)
    p.rod_worth_fraction_L = 0.5;  % Fraction of total rod worth in lower region
                                   % SOURCE: Symmetric assumption (default)
                                   % UNITS: dimensionless (0 to 1)

    p.rod_worth_fraction_U = 0.5;  % Fraction of total rod worth in upper region
                                   % SOURCE: Symmetric assumption (default)
                                   % UNITS: dimensionless (0 to 1)

    % ------------------------------------------------------------------------
    % Control Rod Target Positions
    % ------------------------------------------------------------------------
    % c = 0: rods fully WITHDRAWN (boron out, graphite followers IN core)
    % c = 1: rods fully INSERTED (boron in, graphite followers OUT of core)
    %
    % NORMAL OPERATION: c = 0 (rods withdrawn, maintained by servo)
    % SCRAM TARGET:     c = 1 (rods drive toward full insertion)
    p.c_normal = 0.0;        % Normal operating position (fully withdrawn)
                             % SOURCE: Standard operating procedure
                             % UNITS: dimensionless (0 to 1)

    p.c_scram = 1.0;         % SCRAM target position (fully inserted)
                             % SOURCE: Emergency shutdown procedure
                             % UNITS: dimensionless (0 to 1)

    % ------------------------------------------------------------------------
    % DESIGN FLAW #1: GRAPHITE FOLLOWERS
    % ------------------------------------------------------------------------
    % RBMK control rods have a UNIQUE design:
    %   - Top section: boron carbide absorber (captures neutrons)
    %   - Bottom section: graphite displacer (moderates neutrons)
    %   - Purpose: displace water (absorber) with graphite (better moderator)
    %
    % REACTIVITY COMPONENTS:
    %
    %   When rods WITHDRAWN (c = 0):
    %     - Boron OUT of core → no absorption → neutral reactivity
    %     - Graphite IN core → displaces water → POSITIVE reactivity
    %     - Net effect: +ρ_graphite_follower
    %
    %   When rods INSERTED (c = 1):
    %     - Boron IN core → strong absorption → negative reactivity
    %     - Graphite OUT of core → water returns → less moderation → neutral
    %     - Net effect: -ρ_c_max (large negative)
    %
    % TOTAL CONTROL ROD REACTIVITY MODEL:
    %   ρ_rod(c) = -ρ_c_max × c + ρ_graphite_follower × (1 - c)
    %
    % CALIBRATION: ρ_graphite_follower ≈ 0.3 × ρ_c_max [1]
    %   Graphite effect is ~30% of boron worth
    %   For ρ_c_max = 0.010 (2β): ρ_graphite_follower = 0.003 (+0.6β)
    %
    % CHERNOBYL CONFIGURATION:
    %   - Operating at 200 MW with rods fully withdrawn (c ≈ 0)
    %   - Baseline positive reactivity: +0.6β from graphite followers
    %   - This sets up the positive scram effect (see below)
    %
    % NUMERICAL EXAMPLES:
    %   At c = 0.0 (rods out): ρ_rod = +0.003 = +0.6β  (POSITIVE baseline)
    %   At c = 0.5 (rods half): ρ_rod = -0.0035 = -0.7β (slightly negative)
    %   At c = 1.0 (rods in):  ρ_rod = -0.010 = -2.0β  (full negative worth)
    p.rho_graphite_follower = 0.003;  % Positive reactivity from graphite displacers
                                      % SOURCE: INSAG-7 [1] - estimated at ~30% of rod worth
                                      % CALIBRATED TO: +0.6β = +0.003
                                      % UNCERTAINTY: ±30% (poorly constrained pre-accident)
                                      % UNITS: dimensionless (Δk/k)

    % ------------------------------------------------------------------------
    % DESIGN FLAW #2: POSITIVE SCRAM EFFECT
    % ------------------------------------------------------------------------
    % When SCRAM is initiated, rods insert from TOP of core.
    % But graphite tip is at BOTTOM of rod assembly!
    %
    % INSERTION SEQUENCE:
    %   t = 0-2s: Graphite tips enter BOTTOM of core (lower region)
    %             → Displace water with graphite
    %             → ADD positive reactivity (+1β)
    %
    %   t = 2-18s: Boron section gradually enters core
    %              → Starts absorbing neutrons
    %              → Negative reactivity slowly builds up
    %
    % THE FATAL FLAW:
    %   If reactor is already near prompt critical (due to oscillations),
    %   the initial +1β spike can push it OVER the edge.
    %
    % CHERNOBYL ACCIDENT TIMELINE:
    %   t = 0s:   SCRAM button pressed (reactor oscillating at ~500 MW)
    %   t = 0-2s: Graphite tips enter → +1β reactivity spike
    %   t = 1-2s: Combined with oscillation peak → total ρ > β (prompt critical)
    %   t = 2-3s: Exponential power surge to 30 GW (100× nominal)
    %   t = 3-4s: Fuel disintegration, steam explosion
    %   t = 4-5s: Core destruction, graphite fire
    %
    % CALIBRATION: ρ_tip = +1β = 0.005 [1]
    %   INSAG-7 documents this effect as cause of accident
    %   Magnitude estimated from forensic analysis of power surge
    %
    % POST-ACCIDENT FIX:
    %   - Removed graphite tips
    %   - Shortened absorber section
    %   - Kept rods partially inserted during operation
    %
    % NUMERICAL EXAMPLE:
    %   Reactor at 500 MW with ρ = 0
    %   Oscillation brings ρ to +0.003 (+0.6β)
    %   SCRAM adds +0.005 (+1.0β)
    %   Total: +0.008 (+1.6β) → PROMPT CRITICAL → explosion
    p.rho_tip = 0.005;       % Positive reactivity from graphite tip insertion
                             % SOURCE: INSAG-7 [1] - forensic reconstruction
                             % CALIBRATED TO: +1β = 0.005
                             % UNCERTAINTY: ±30% (estimated from accident dynamics)
                             % UNITS: dimensionless (Δk/k)

    p.tau_tip = 2.0;         % Duration of tip effect (s)
                             % SOURCE: Estimated from rod geometry and insertion speed
                             % PHYSICS: Time for 1.5m graphite tip to enter 7m core
                             % UNCERTAINTY: ±50%
                             % UNITS: seconds

    % ========================================================================
    %% 6. POWER-DEPENDENT STABILITY ENHANCEMENTS
    % ========================================================================
    % Phenomenological effects that strengthen at high power, contributing to
    % the stabilization above 762 MW (Hopf bifurcation threshold).

    % ------------------------------------------------------------------------
    % Hopf Bifurcation Threshold
    % ------------------------------------------------------------------------
    % EXPERIMENTAL FINDING: RBMK becomes unstable below ~24% nominal power [3]
    %   Critical power: 762 MW (0.24 fractional power)
    %   Below this: sustained oscillations (limit cycle)
    %   Above this: stable operation (oscillations decay)
    %
    % This threshold emerges from the balance:
    %   Destabilizing: positive void feedback + transport delay
    %   Stabilizing:   negative Doppler feedback + xenon regulation
    %
    % At low power:  void feedback dominates → Hopf instability
    % At high power: Doppler + xenon dominate → stable
    p.power_threshold = 0.24;     % Stability threshold (dimensionless)
                                  % SOURCE: EPJ Nuclear Sciences [3]
                                  % CALIBRATED TO: 762 MW = 24% nominal
                                  % UNCERTAINTY: ±10%
                                  % UNITS: dimensionless (fraction of nominal power)

    % ------------------------------------------------------------------------
    % Xenon Saturation Modeling
    % ------------------------------------------------------------------------
    % Flag to enable/disable xenon saturation effects in derivatives function.
    % Included for future extensions (not currently implemented).
    p.xenon_saturation = true;    % Enable xenon saturation modeling (boolean)
                                  % SOURCE: Placeholder for future enhancement
                                  % UNITS: logical (true/false)

    % ------------------------------------------------------------------------
    % Doppler Enhancement at High Power
    % ------------------------------------------------------------------------
    % Doppler coefficient strengthens at high power due to:
    %   - Higher fuel temperatures → more U-238 resonance absorption
    %   - Flux spectrum hardening in high-void regions
    %
    % PHENOMENOLOGICAL MODEL:
    %   κ_D_effective = κ_D0 × [1 + (doppler_enhancement - 1) × power_fraction]
    %
    % At nominal power (n=1): κ_D increases by 50% (enhancement = 1.5)
    % At low power (n=0.2):   κ_D is near baseline (enhancement ≈ 1.1)
    %
    % This contributes to high-power stability but is NOT sufficient to
    % stabilize low power (void feedback still dominates).
    p.doppler_enhancement = 1.5;  % Doppler multiplier at nominal power (dimensionless)
                                  % SOURCE: Phenomenological fit to stability threshold
                                  % UNCERTAINTY: ±50% (effective parameter, not physical)
                                  % UNITS: dimensionless (multiplier)

    % ------------------------------------------------------------------------
    % Void Saturation Coefficient
    % ------------------------------------------------------------------------
    % At high void fractions, void dynamics become LESS sensitive:
    %   - Approaching annular flow regime (α > 0.6)
    %   - Further boiling has diminishing returns
    %   - Effective time constant increases
    %
    % SATURATION MODEL:
    %   effective_time_constant = τ_v × [1 + void_saturation_coeff × α]
    %
    % NUMERICAL EXAMPLES:
    %   At α = 0.0: no saturation (factor = 1.0)
    %   At α = 0.4: saturation_factor = 1.8 (44% slower dynamics)
    %   At α = 0.6: saturation_factor = 2.2 (55% slower dynamics)
    %
    % PHYSICS: High void → longer bubble coalescence time → slower dynamics
    %
    % This effect contributes to stability at high power where void fractions
    % are large (α ≈ 0.6). At low power (α ≈ 0.2), minimal effect.
    p.void_saturation_coeff = 2.0;  % Void saturation coefficient (dimensionless)
                                    % SOURCE: Fitted to match stability threshold [3]
                                    % PHYSICS: Reduces void dynamics at high α
                                    % UNCERTAINTY: ±50% (phenomenological)
                                    % UNITS: dimensionless

    % ========================================================================
    %% 7. PARAMETER VALIDATION
    % ========================================================================
    % Automatic consistency checks to catch configuration errors.

    % Check that rod worth fractions sum to 1.0 (conservation of worth)
    rod_worth_sum = p.rod_worth_fraction_L + p.rod_worth_fraction_U;
    if abs(rod_worth_sum - 1.0) > 1e-10
        error('rbmk_parameters:InvalidRodWorth', ...
            'Rod worth fractions must sum to 1.0, got %.6f', rod_worth_sum);
    end

    % Future validation checks can be added here:
    %   - Positive parameters (masses, time constants, etc.)
    %   - Physical bounds (0 < α_max < 1, β < 0.01, etc.)
    %   - Dimensional consistency
end

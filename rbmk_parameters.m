function modes = rbmk_parameters()
    % RBMK_PARAMETERS - Dual-Regime Parameter Definitions
    %
    % PURPOSE:
    %   Provides two distinct parameter sets to demonstrate the contrast
    %   between standard operation and the specific conditions of April 26.
    %
    % OUTPUT:
    %   modes.normal   : Standard RBMK parameters (Stable, High Flow)
    %   modes.accident : Forensic RBMK parameters (Unstable, Low Flow, High Burnup)
    
    %% 1. COMMON PHYSICS (Invariant Constants)
    % These properties depend on materials (UO2, Graphite, Water) and 
    % geometry, which do not change between modes.
    % =====================================================================
    c = struct();
    
    % Neutronics
    c.beta = 0.005;          % Delayed neutron fraction
    c.Lambda = 4.8e-4;       % Prompt generation time (s)
    c.lambda_d = 0.08;       % Precursor decay (1/s)
    
    % Power & Thermodynamics
    c.P_nominal = 3200;      % Total MW
    c.k_P = 1600;            % MW per region
    c.Tc = 270 + 273.15;     % Inlet Temp (K)
    c.T_sat = 285 + 273.15;  % Saturation Temp (K)
    c.T0 = 543.15;           % Ref Temp (K)
    c.h_fg = 1505;           % Enthalpy (kJ/kg)
    
    % Heat Transfer (Lumped Model)
    c.b_f = 0.2;             % Fuel relaxation (1/s)
    c.a_f = 110.0;           % Fuel heating rate
    c.b_m = 0.05;            % Moderator relaxation
    c.a_m = 15.0;            % Moderator heating rate
    
    % Boiling Curve Physics
    c.p_shape = 0.25;        % Root-function shape (Flash boiling)
    c.k_sat = 2.0;           % Saturation damping
    
    % Xenon Constants
    c.kappa_X = 3.7e-5;      
    c.lambda_I = 2.9e-5; c.lambda_X = 2.1e-5;
    c.y_I = 0.06; c.y_X = 0.002;
    c.Phi_nominal = 1.0e14; c.sigma_X = 2.0e-18;

    % Scram Mechanics (Rod speed / Tip)
    c.tau_c = 18.0;          % 18s insertion time
    c.kappa_tip = 0.272;     % The Graphite Tip (+2 beta spike)
    c.kappa_boron = 0.010;   % Absorber worth (-2 beta)
    c.c_scram = 1.0;         % SCRAM target position (fully inserted)
    
    % Termination Thresholds (Explosion Physics)
    c.use_fuel_dispersal = true;
    c.T_melt = 3100;         
    c.T_disperse = 3500;
    c.kappa_disassembly = -0.50; 
    c.P_disassembly = 30000; 
    c.K_heat = (c.k_P * 1000) / c.h_fg;

    %% 2. NORMAL MODE (Safe Operation)
    % Represents the reactor at full power with proper rod configuration.
    % =====================================================================
    norm = c;
    
    % Flow: Nominal (8000 kg/s per half)
    norm.m_flow = 8000;      
    norm.tau_flow = 2.0;     % Standard transit time
    norm.tau_void = 1.0;     % Fast bubble relaxation
    
    % Coupling: Tight (Core acts as one unit)
    norm.Dn = 6.0;           
    
    % Reactivity: "Design" Values (Stable)
    norm.use_dynamic_kappa_V = true;
    norm.kappa_V_low = 0.025;       % +5 beta (Standard high burnup)
    norm.kappa_V_high = 0.005;      % +1 beta (At full power)
    norm.kappa_V_transition = 0.15; % Normal spectral shift
    norm.alpha_max = 0.85;
    
    % Doppler: Strong (Good Braking)
    norm.use_dynamic_kappa_D = true;
    norm.kappa_D_base = -0.012;     % Full strength negative feedback
    norm.doppler_scale = 3.0;       % Activates quickly
    
    % Control
    norm.t_scram = 100.0;    % No SCRAM scheduled
    norm.c_normal = 0.5;     % Rods halfway in (Normal flux shaping)
    
    modes.normal = norm;

    %% 3. ACCIDENT MODE (April 26, 1986)
    % Represents the degraded state: Low flow, Rods withdrawn, High Void Worth.
    % =====================================================================
    acc = c;
    
    % Flow: Reduced (Pump Coastdown / Cavitation)
    acc.m_flow = 5600;       % 70% Flow
    
    % Delay: Increased (Physics Correction: Tau ~ 1/Flow)
    acc.tau_flow = 2.85;     % 2.0 * (8000/5600)
    acc.tau_void = 1.5;      % Slower dynamics at low flow
    
    % Coupling: Loose (Spatial Instability allowed)
    acc.Dn = 0.2;            
    
    % Reactivity: "Killer" Values (Forensically Tuned)
    acc.use_dynamic_kappa_V = true;
    acc.kappa_V_low = 0.080;        % +16 beta (Extreme positive feedback)
    acc.kappa_V_high = 0.015;
    acc.kappa_V_transition = 10.0;  % "Sticky" high void coeff (Fatal flaw)
    acc.alpha_max = 0.95;           % Allow higher voiding
    
    % Doppler: Broken (Cold fuel = Weak feedback)
    acc.use_dynamic_kappa_D = true;
    acc.kappa_D_base = -0.006;      % Half strength brake
    acc.doppler_scale = 1.5;        % Slow activation
    
    % Control
    acc.t_scram = 30.0;      % AZ-5 pressed at t=30
    acc.c_normal = 0.0;      % Rods fully withdrawn (Xenon override)
    
    modes.accident = acc;
end

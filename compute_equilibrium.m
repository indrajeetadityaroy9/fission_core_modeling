function [y_eq, p_corrected, diagnostics] = compute_equilibrium(target_power_fraction, p, options)
    % COMPUTE_EQUILIBRIUM_ROBUST - Robust Initializer for RBMK Model
    %
    % DESCRIPTION:
    %   Calculates the steady-state vector y_eq (16x1) for a given power level.
    %   Uses Constrained Minimization (fminbnd) to solve the criticality equation,
    %   ensuring robustness against the non-linear "Tip Effect" and "Xenon Pit".
    %
    % ROBUSTNESS UPGRADES:
    %   1. SOLVER: Replaced fzero with fminbnd. This prevents the solver from
    %      wandering outside physical bounds [0, 1] or crashing if no root exists.
    %   2. PHYSICS: Automatically computes 'rho_0' (Excess Reactivity) if the
    %      control rods alone are insufficient to balance the core.
    %   3. UNITS: Strictly adheres to Kelvin temperatures and Flux scaling.
    %
    % INPUTS:
    %   target_power_fraction - Normalized power (0.0 to 1.0)
    %   p                     - Parameter struct from rbmk_parameters()
    %   options               - (Optional) 'normal' or 'accident'
    %
    % OUTPUTS:
    %   y_eq        - 16x1 State Vector [L; U] ready for DDE solver
    %   p_corrected - Parameter struct containing the calculated 'rho_0'
    %   diagnostics - Struct with reactivity breakdown
    %
    % AUTHOR: RBMK Forensics Model (Optimized Version)

    %% 1. INITIALIZATION & SETUP
    if nargin < 3, options = 'normal'; end
    accident_mode = strcmpi(options, 'accident');

    % Check if user passed the container struct instead of a specific mode
    % rbmk_parameters() returns modes.normal and modes.accident
    if isfield(p, 'normal') && isfield(p, 'accident')
        % User passed the container - extract the appropriate mode
        if accident_mode
            p = p.accident;
            warning('compute_equilibrium:ContainerDetected', ...
                'Received parameter container. Auto-selecting p.accident based on options.');
        else
            p = p.normal;
            warning('compute_equilibrium:ContainerDetected', ...
                'Received parameter container. Auto-selecting p.normal based on options.');
        end
    end

    p_corrected = p;
    % Ensure rho_0 exists (will be overwritten if needed)
    if ~isfield(p_corrected, 'rho_0'), p_corrected.rho_0 = 0.0; end

    % Set Flux for Xenon Burnout (Must match p.Phi_nominal logic)
    if isfield(p, 'Phi_nominal')
        Phi = p.Phi_nominal;
    else
        Phi = 1.0e14;
    end

    %% 1.5. PARAMETER INTEGRITY CHECK
    % Verify that the selected 'p' struct has the required physics fields.
    % This catches errors where the user passes an empty or malformed struct.

    required_fields = {'beta', 'Lambda', 'm_flow', 'tau_flow', 'Dn', ...
                       'K_heat', 'kappa_X', 'lambda_I', 'lambda_X', ...
                       'c_scram', 'c_normal', 't_scram', 'tau_c'};

    for i = 1:length(required_fields)
        if ~isfield(p, required_fields{i})
            error('compute_equilibrium:MissingField', ...
                'Parameter struct is missing required field: %s', required_fields{i});
        end
    end

    % Verify Flow/Delay Consistency (Physics Sanity Check)
    % Nominal: 8000 kg/s -> 2.0s, so tau_flow ~ 16000 / m_flow
    estimated_tau = 16000 / p.m_flow;
    if abs(p.tau_flow - estimated_tau) > 0.5
        warning('compute_equilibrium:PhysicsMismatch', ...
            'tau_flow (%.2fs) does not match flow scaling (approx %.2fs). Check parameters.', ...
            p.tau_flow, estimated_tau);
    end

    %% 2. CALCULATE ALGEBRAIC STATE VARIABLES
    % Solve differential equations at steady state (dy/dt = 0)
    
    n = target_power_fraction;

    % A. Precursors: C = (beta / (Lambda * lambda_d)) * n
    C = (p.beta / (p.Lambda * p.lambda_d)) * n;

    % B. Temperatures (Kelvin): T = Tc + (a/b)*n
    Tf = p.Tc + (p.a_f / p.b_f) * n; 
    Tm = p.Tc + (p.a_m / p.b_m) * n;

    % C. Iodine & Xenon (with Flux Burnout)
    I = (p.y_I / p.lambda_I) * n;
    
    production = (p.y_X * n) + (p.lambda_I * I);
    burnout_rate = p.sigma_X * Phi * n; % Rate in 1/s
    X = production / (p.lambda_X + burnout_rate);

    % D. Thermal-Hydraulics (Boiling Curve)
    x_quality = (p.K_heat * n) / p.m_flow;
    x_quality = max(0, x_quality + 1e-9); % Prevent numerical zero
    alpha = p.alpha_max * (x_quality^p.p_shape) / (1 + x_quality^p.p_shape);

    %% 3. CALCULATE INHERENT REACTIVITY
    % These feedbacks are fixed by the physics state (n, T, alpha, X)
    
    % Void Feedback (Dynamic check)
    if isfield(p, 'use_dynamic_kappa_V') && p.use_dynamic_kappa_V
        n_safe = max(n, 0.01);
        kappa_V_eq = p.kappa_V_low - (p.kappa_V_low - p.kappa_V_high) * ...
                     (1 - exp(-n_safe / p.kappa_V_transition));
    else
        kappa_V_eq = p.kappa_V;
    end
    rho_void = kappa_V_eq * alpha;

    % Doppler Feedback (Dynamic check, using Kelvin T0)
    if isfield(p, 'use_dynamic_kappa_D') && p.use_dynamic_kappa_D
        n_safe = max(n, 0.01);
        kappa_D_eq = p.kappa_D_base * (1 - exp(-n_safe * p.doppler_scale));
    else
        kappa_D_eq = p.kappa_D;
    end
    rho_dop = kappa_D_eq * (sqrt(Tf) - sqrt(p.T0));

    % Xenon Feedback
    rho_xen = -p.kappa_X * X;

    % Total Inherent Reactivity (The "Load" the rods must carry)
    rho_inherent = rho_void + rho_dop + rho_xen;

    %% 4. CRITICALITY SEARCH (ROBUST OPTIMIZATION)
    % Goal: Find c and rho_0 such that: rho_inherent + rho_rod(c) + rho_0 = 0
    
    % Rod Reactivity Function (Lower Core Physics includes Tip Effect)
    rod_func = @(c) (p.kappa_tip .* c .* exp(-10.*c) - p.kappa_boron .* c);

    if accident_mode
        % --- ACCIDENT SCENARIO ---
        % Force rods fully withdrawn (c=0) to simulate the "Trap".
        % Calculate exactly how much excess reactivity is needed to be critical.
        c_eq = 0.0;
        
        % At c=0, rho_rod is 0.
        % The missing reactivity must come from rho_0.
        rho_missing = -(rho_inherent + rod_func(c_eq));
        
        p_corrected.rho_0 = rho_missing;
        status = 'ACCIDENT MODE (c=0 forced)';
        
    else
        % --- NORMAL OPERATION ---
        % Try to find a rod position c in [0, 1] that balances the core.
        % We minimize the squared error of reactivity.
        
        rho_target = -rho_inherent;
        
        % Objective: Minimize |rho_rod(c) - rho_target|
        objective = @(c) abs(rod_func(c) - rho_target);
        
        % fminbnd is robust: it respects bounds and handles non-roots.
        opts = optimset('TolX', 1e-8, 'Display', 'off');
        [c_eq, error_val] = fminbnd(objective, 0.0, 1.0, opts);
        
        % Calculate the residual (what the rods couldn't handle)
        rho_rod_actual = rod_func(c_eq);
        residual = rho_inherent + rho_rod_actual;
        
        % If error is small, rods balanced it. If large, add rho_0.
        if error_val < 1e-7
            p_corrected.rho_0 = 0.0;
            status = 'Converged (Rods only)';
        else
            % Rods hit a limit (0 or 1), use rho_0 to make up the rest.
            p_corrected.rho_0 = -residual;
            status = sprintf('Limit Reached (c=%.2f, rho_0 added)', c_eq);
        end
    end

    %% 5. ASSEMBLE OUTPUT VECTORS
    % Initialize both regions symmetrically
    % State Vector: [n, C, alpha, Tf, Tm, I, X, c]
    y_L = [n; C; alpha; Tf; Tm; I; X; c_eq];
    y_U = [n; C; alpha; Tf; Tm; I; X; c_eq];
    
    y_eq = [y_L; y_U];

    %% 6. DIAGNOSTICS
    diagnostics.n = n;
    diagnostics.c_eq = c_eq;
    diagnostics.rho_0 = p_corrected.rho_0;
    diagnostics.rho_void = rho_void;
    diagnostics.rho_dop = rho_dop;
    diagnostics.rho_xen = rho_xen;
    diagnostics.rho_rod = rod_func(c_eq);
    diagnostics.rho_total = rho_inherent + rod_func(c_eq) + p_corrected.rho_0;
    diagnostics.status = status;
    
    % Verification Printout
    fprintf('COMPUTE_EQUILIBRIUM [Power=%.2f]: %s\n', n, status);
    fprintf('  Reactivity Breakdown:\n');
    fprintf('    Void:    %+.5f\n', rho_void);
    fprintf('    Doppler: %+.5f\n', rho_dop);
    fprintf('    Xenon:   %+.5f\n', rho_xen);
    fprintf('    Rod(c):  %+.5f (c=%.4f)\n', rod_func(c_eq), c_eq);
    fprintf('    Excess:  %+.5f (rho_0)\n', p_corrected.rho_0);
    fprintf('    TOTAL:   %+.1e (Should be ~0)\n', diagnostics.rho_total);
end

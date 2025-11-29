function [y_eq, p_corrected, diagnostics] = compute_equilibrium(target_power_fraction, p)
    % COMPUTE_EQUILIBRIUM - Robust steady-state solver with Criticality Search
    %
    % DESCRIPTION:
    %   Performs a "Critical Reactivity Search" to initialize the RBMK model.
    %   If the physics parameters (Doppler/Xenon) overwhelm the control rods,
    %   this function automatically calculates a 'bias' reactivity (rho_0)
    %   to force criticality, simulating the core's excess fuel reactivity.
    %
    %   This implements the fundamental reactor physics principle that fresh
    %   fuel must have EXCESS REACTIVITY to overcome:
    %     - Doppler broadening (negative feedback from fuel heating)
    %     - Xenon-135 poisoning (negative feedback from fission products)
    %     - Control rod absorption (negative reactivity for power control)
    %
    % USAGE:
    %   [y0, p_new, diag] = compute_equilibrium(0.5, p);  % 50% power
    %   % IMPORTANT: Pass p_new to rbmk_dynamics, not the original p!
    %
    % INPUTS:
    %   target_power_fraction - Normalized power (0.0 to 1.0)
    %                           n=1.0 corresponds to p.P_nominal per region
    %   p                     - Parameter structure from rbmk_parameters()
    %
    % OUTPUTS:
    %   y_eq        - 16x1 State Vector ready for DDE integration
    %   p_corrected - Parameter struct with added field 'rho_0' (Excess Reactivity)
    %                 YOU MUST USE THIS for simulation, not the original p!
    %   diagnostics - Structure containing reactivity breakdown for verification
    %
    % PHYSICS:
    %   At steady state, total reactivity must be zero:
    %     rho_total = rho_void + rho_doppler + rho_xenon + rho_rod + rho_0 = 0
    %
    %   Where rho_0 represents the excess reactivity from fuel enrichment.
    %   This function solves for rho_0 when the control rod range is insufficient.
    %
    % Author: RBMK Forensics Model

    %% 1. INITIALIZE & COPY PARAMETERS
    p_corrected = p;
    if ~isfield(p_corrected, 'rho_0')
        p_corrected.rho_0 = 0.0; % Initialize if not present
    end

    % Flux Scaling for Xenon burnout
    if isfield(p, 'Phi_nominal')
        Phi = p.Phi_nominal;
    else
        Phi = 1.0e14;  % Default RBMK thermal flux (n/cm^2/s)
    end

    %% 2. CALCULATE STATE VARIABLES (Algebraic Equilibrium)
    n = target_power_fraction;

    % ------------------------------------------------------------------------
    % A. NEUTRON KINETICS: Precursor Equilibrium
    % ------------------------------------------------------------------------
    % From dC/dt = (beta/Lambda)*n - lambda_d*C = 0
    % Solving: C_eq = (beta / (Lambda * lambda_d)) * n
    %
    % For RBMK: ratio = 0.005 / (4.8e-4 * 0.08) = 130.2
    % CRITICAL: C must be ~130x larger than n, or power will crash!
    C = (p.beta / (p.Lambda * p.lambda_d)) * n;

    % ------------------------------------------------------------------------
    % B. THERMAL EQUILIBRIUM
    % ------------------------------------------------------------------------
    % From dT/dt = a*n - b*(T - Tc) = 0
    % Solving: T_eq = Tc + (a/b)*n
    Tf = p.Tc + (p.a_f / p.b_f) * n;  % Fuel temperature (K)
    Tm = p.Tc + (p.a_m / p.b_m) * n;  % Moderator temperature (K)

    % ------------------------------------------------------------------------
    % C. IODINE & XENON EQUILIBRIUM
    % ------------------------------------------------------------------------
    % Iodine: dI/dt = y_I*n - lambda_I*I = 0
    I = (p.y_I / p.lambda_I) * n;

    % Xenon: dX/dt = y_X*n + lambda_I*I - (lambda_X + sigma_X*Phi*n)*X = 0
    % Production = Direct yield + Iodine decay
    % Loss = Radioactive decay + Neutron burnout (flux-dependent!)
    production = (p.y_X * n) + (p.lambda_I * I);
    burnout_rate = p.sigma_X * Phi * n;  % [cm^2] * [1/cm^2/s] * [1] = [1/s]
    X = production / (p.lambda_X + burnout_rate);

    % ------------------------------------------------------------------------
    % D. THERMAL-HYDRAULIC EQUILIBRIUM (Void Fraction)
    % ------------------------------------------------------------------------
    % Steam quality: x = K_heat * n / m_flow
    x = (p.K_heat * n) / p.m_flow;
    x = max(0, x + 1e-9);  % Small epsilon for numerical stability

    % Void fraction from boiling curve (equilibrium, no saturation factor)
    alpha = p.alpha_max * (x^p.p_shape) / (1 + x^p.p_shape);

    %% 3. CALCULATE INHERENT REACTIVITY FEEDBACKS
    rho_void = p.kappa_V * alpha;                           % Positive (RBMK flaw)
    rho_dop  = p.kappa_D * (sqrt(Tf) - sqrt(p.T0));        % Negative (stabilizing)
    rho_xen  = -p.kappa_X * X;                              % Negative (poison)

    % Net inherent reactivity (without rods or excess fuel)
    rho_inherent = rho_void + rho_dop + rho_xen;

    %% 4. CRITICALITY SEARCH (The Critical Step)
    % We need: rho_inherent + rho_rod(c) + rho_0 = 0
    %
    % Strategy:
    %   1. First assume rho_0 = 0. Can control rods balance the core?
    %   2. If yes, solve for c_eq using fzero.
    %   3. If no (rods insufficient), calculate required rho_0.

    % Define Rod Reactivity Function (with tip effect)
    % Use element-wise operations (.*) for vectorized evaluation
    rod_func = @(c) (p.kappa_tip .* c .* exp(-10.*c) - p.kappa_boron .* c);

    % We need rho_rod = -rho_inherent to achieve criticality
    rho_target = -rho_inherent;

    % Determine the achievable range of rod reactivity
    c_vals = linspace(0, 1, 100);
    rho_vals = rod_func(c_vals);
    min_rod = min(rho_vals);  % Most negative (rods fully in)
    max_rod = max(rho_vals);  % Most positive (tip effect peak)

    if rho_target >= min_rod && rho_target <= max_rod
        % CASE A: Control rods can achieve criticality
        % Solve rod_func(c) = rho_target
        %
        % CHALLENGE: rod_func is NON-MONOTONIC (peaks at c≈0.1)
        % fzero needs an interval where the function changes sign.
        %
        % Strategy: Find the peak, then search appropriate interval
        [~, peak_idx] = max(rho_vals);
        c_peak = c_vals(peak_idx);
        rho_peak = max_rod;

        problem = @(c) rod_func(c) - rho_target;

        try
            if rho_target >= 0 && rho_target <= rho_peak
                % Target is positive: solution is on the rising side (0 to peak)
                % Check if there's a sign change
                if problem(0) * problem(c_peak) <= 0
                    c_eq = fzero(problem, [0, c_peak]);
                else
                    % No sign change, target might equal peak or endpoint
                    c_eq = fzero(problem, c_peak);  % Start from peak
                end
            else
                % Target is negative: solution is on the falling side (peak to 1)
                if problem(c_peak) * problem(1) <= 0
                    c_eq = fzero(problem, [c_peak, 1]);
                else
                    c_eq = fzero(problem, 0.5);  % Fallback starting point
                end
            end
            p_corrected.rho_0 = 0.0;
            status = 'Converged (Rods only)';
        catch
            % Fallback: set c to achieve closest match, add rho_0 for balance
            if rho_target > 0
                c_eq = c_peak;  % Use peak for max positive reactivity
            else
                c_eq = 1.0;     % Use full insertion for max negative
            end
            p_corrected.rho_0 = -(rho_inherent + rod_func(c_eq));
            status = sprintf('Fallback (c=%.2f, rho_0=%+.4f)', c_eq, p_corrected.rho_0);
        end
    else
        % CASE B: Rods cannot balance the core (need excess reactivity)
        % This is normal! Real reactors have excess fuel reactivity.
        %
        % Set rods to "normal" position (withdrawn, c=0)
        c_eq = 0.0;

        % Calculate the missing reactivity (excess fuel reactivity)
        % Criticality: 0 = rho_inherent + rho_rod(0) + rho_0
        rho_rod_actual = rod_func(c_eq);
        p_corrected.rho_0 = -(rho_inherent + rho_rod_actual);

        status = sprintf('Converged (rho_0 = %+.4f = %+.2f beta)', ...
            p_corrected.rho_0, p_corrected.rho_0/p.beta);
    end

    %% 5. ASSEMBLE STATE VECTOR
    % Both regions start in symmetric equilibrium
    y_L = [n; C; alpha; Tf; Tm; I; X; c_eq];
    y_U = [n; C; alpha; Tf; Tm; I; X; c_eq];
    y_eq = [y_L; y_U];

    %% 6. DIAGNOSTICS OUTPUT
    diagnostics.rho_void = rho_void;
    diagnostics.rho_dop = rho_dop;
    diagnostics.rho_xen = rho_xen;
    diagnostics.rho_rod = rod_func(c_eq);
    diagnostics.rho_0 = p_corrected.rho_0;
    diagnostics.rho_total = rho_void + rho_dop + rho_xen + rod_func(c_eq) + p_corrected.rho_0;
    diagnostics.c_eq = c_eq;
    diagnostics.void_fraction = alpha;
    diagnostics.fuel_temp = Tf;
    diagnostics.xenon = X;
    diagnostics.precursor_ratio = C / n;
    diagnostics.status = status;

    fprintf('Equilibrium [n=%.3f]: %s\n', n, status);
    fprintf('  Reactivity: void=%+.4f, dop=%+.4f, xen=%+.4f, rod=%+.4f, rho_0=%+.4f\n', ...
        rho_void, rho_dop, rho_xen, rod_func(c_eq), p_corrected.rho_0);
    fprintf('  Total rho = %.2e (should be ~0)\n', diagnostics.rho_total);
end

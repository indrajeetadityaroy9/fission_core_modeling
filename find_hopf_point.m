function [P_boundary, omega_boundary, info] = find_hopf_point(power_levels, eigenvalues_DDE, varargin)
    % ANALYZE_STABILITY_BOUNDARY - Detect stability boundary from eigenvalue data
    %
    % Identifies the stability boundary by finding where eigenvalues cross the
    % imaginary axis as power varies. This detects the point where the system
    % transitions from stable to unstable (or vice versa).
    %
    % THEORY:
    %   Linear stability boundary occurs when:
    %     - Dominant eigenvalue λ = α(P) ± iω(P)
    %     - Real part α crosses zero: α < 0 (stable) → α > 0 (unstable)
    %     - At boundary: α = 0
    %
    % BIFURCATION TYPES:
    %   - Hopf bifurcation: Complex pair crosses axis (ω ≠ 0) → oscillatory instability
    %   - Saddle-node/divergent: Real eigenvalue crosses (ω = 0) → exponential instability
    %
    % FOR RBMK REACTOR:
    %   This function detects the DIVERGENT INSTABILITY BOUNDARY at ~1475 MW.
    %   Below this power, the system is unstable (eigenvalue real part > 0).
    %   The instability is NOT oscillatory (Hopf) but EXPONENTIAL (divergent).
    %
    %   Physical interpretation:
    %     - Power > 1475 MW: Stable operation (α < 0, perturbations decay)
    %     - Power < 1475 MW: Divergently unstable (α > 0, exponential growth)
    %     - Accident at 200 MW: Deep in unstable region, exponential power surge
    %
    % INPUTS:
    %   power_levels    - Vector of power levels (MW) [n_points × 1]
    %   eigenvalues_DDE - Cell array of eigenvalue vectors {n_points × 1}
    %                     Each cell contains eigenvalues at that power level
    %
    % OPTIONAL PARAMETERS (Name-Value pairs):
    %   'MinFrequency'  - Minimum oscillation frequency to consider (default: 0.01 rad/s)
    %   'Verbose'       - Display detailed diagnostic info (default: false)
    %
    % OUTPUTS:
    %   P_boundary     - Power level at stability boundary (MW)
    %   omega_boundary - Oscillation frequency at boundary (rad/s)
    %                    = 0 for divergent instability (RBMK case)
    %                    ≠ 0 for Hopf bifurcation (oscillatory instability)
    %   info           - Structure with:
    %                    .eigenvalue_trace - Real parts of dominant mode vs power
    %                    .frequency_trace - Imaginary parts vs power
    %                    .crossing_index - Index where zero crossing occurs
    %                    .interpolation_quality - Quality metric for interpolation
    %                    .detection_method - Method used to find boundary
    %                    .instability_type - 'divergent' or 'hopf'
    %
    % EXAMPLE:
    %   % After running hopf_analysis
    %   p = rbmk_parameters();
    %   [powers, amps, periods, outcomes, eig_J0, eig_DDE, margin_J0, margin_DDE] = hopf_analysis(p);
    %   [P_boundary, omega, info] = find_hopf_point(powers, eig_DDE);
    %
    %   fprintf('Stability boundary at %.0f MW\n', P_boundary);
    %   if omega < 0.01
    %       fprintf('Instability type: Divergent (exponential growth)\n');
    %   else
    %       fprintf('Instability type: Hopf (oscillatory, period %.1f s)\n', 2*pi/omega);
    %   end
    %
    % See also: hopf_analysis, compute_stability

    %% Parse inputs
    p_input = inputParser();
    addRequired(p_input, 'power_levels', @(x) isnumeric(x) && isvector(x));
    addRequired(p_input, 'eigenvalues_DDE', @(x) iscell(x));
    addParameter(p_input, 'MinFrequency', 0.01, @(x) isnumeric(x) && isscalar(x));
    addParameter(p_input, 'Verbose', false, @islogical);
    parse(p_input, power_levels, eigenvalues_DDE, varargin{:});

    opts = p_input.Results;
    power_levels = power_levels(:);  % Ensure column vector
    n_points = length(power_levels);

    if opts.Verbose
        fprintf('=== STABILITY BOUNDARY ANALYSIS ===\n');
        fprintf('Analyzing %d power levels from %.0f to %.0f MW\n\n', ...
            n_points, min(power_levels), max(power_levels));
    end

    %% Extract dominant mode at each power level
    eigenvalue_trace = zeros(n_points, 1);
    frequency_trace = zeros(n_points, 1);

    for i = 1:n_points
        eigs_i = eigenvalues_DDE{i};
        if isempty(eigs_i)
            eigenvalue_trace(i) = NaN;
            frequency_trace(i) = NaN;
            continue;
        end

        % Find complex eigenvalues with significant imaginary part
        complex_idx = abs(imag(eigs_i)) > opts.MinFrequency;
        complex_eigs = eigs_i(complex_idx);

        if isempty(complex_eigs)
            % No oscillatory modes, use rightmost eigenvalue (divergent case)
            eigenvalue_trace(i) = real(eigs_i(1));
            frequency_trace(i) = 0;
        else
            % Use complex pair with largest real part (Hopf case)
            [~, max_idx] = max(real(complex_eigs));
            dominant_mode = complex_eigs(max_idx);

            eigenvalue_trace(i) = real(dominant_mode);
            frequency_trace(i) = abs(imag(dominant_mode));
        end
    end

    %% Find zero crossing of real part
    valid_idx = ~isnan(eigenvalue_trace);
    valid_powers = power_levels(valid_idx);
    valid_reals = eigenvalue_trace(valid_idx);
    valid_freqs = frequency_trace(valid_idx);

    if sum(valid_idx) < 2
        error('find_hopf_point:InsufficientData', ...
            'Need at least 2 valid eigenvalue points');
    end

    % Find sign changes (going from high power to low power: negative → positive)
    % ENHANCEMENT: Add tolerance-based detection for near-zero crossings
    zero_tolerance = 1e-4;  % Threshold for "essentially zero" margin

    % Method 1: Traditional sign change detection
    sign_changes = find(diff(sign(valid_reals)) ~= 0);

    % Method 2: Tolerance-based crossing (for near-zero margins)
    near_zero_crossings = [];
    for i = 1:length(valid_reals)-1
        % Check if we transition through zero within tolerance
        if abs(valid_reals(i)) < zero_tolerance && abs(valid_reals(i+1)) < zero_tolerance
            % Both points near zero - potential crossing
            if sign(valid_reals(i)) ~= sign(valid_reals(i+1))
                near_zero_crossings = [near_zero_crossings, i];
            end
        elseif abs(valid_reals(i)) < zero_tolerance || abs(valid_reals(i+1)) < zero_tolerance
            % One point very close to zero
            near_zero_crossings = [near_zero_crossings, i];
        end
    end

    % Method 3: Fallback - find point closest to zero
    [min_abs_real, closest_idx] = min(abs(valid_reals));

    % Decide which method to use
    if ~isempty(sign_changes)
        % Traditional sign change found
        crossing_idx = sign_changes(1);
        detection_method = 'sign_change';
    elseif ~isempty(near_zero_crossings)
        % Near-zero crossing found
        crossing_idx = near_zero_crossings(1);
        detection_method = 'near_zero';
        warning('find_hopf_point:NearZeroCrossing', ...
            'No clear sign change, using near-zero tolerance detection');
    elseif min_abs_real < zero_tolerance
        % Fallback: use closest-to-zero point
        if closest_idx == 1
            crossing_idx = 1;
        elseif closest_idx == length(valid_reals)
            crossing_idx = length(valid_reals) - 1;
        else
            % Use interval that straddles the minimum
            if abs(valid_reals(closest_idx-1)) < abs(valid_reals(closest_idx+1))
                crossing_idx = closest_idx - 1;
            else
                crossing_idx = closest_idx;
            end
        end
        detection_method = 'closest_to_zero';
        warning('find_hopf_point:NoZeroCrossing', ...
            'No zero crossing found, using closest-to-zero point (margin = %.6f)', min_abs_real);
    else
        % No crossing found by any method
        warning('find_hopf_point:NoZeroCrossing', ...
            'No zero crossing found in eigenvalue real part (min |margin| = %.6f)', min_abs_real);
        P_boundary = NaN;
        omega_boundary = NaN;
        info.eigenvalue_trace = eigenvalue_trace;
        info.frequency_trace = frequency_trace;
        info.crossing_index = [];
        info.interpolation_quality = NaN;
        info.detection_method = 'failed';
        info.min_abs_margin = min_abs_real;
        info.instability_type = 'unknown';
        return;
    end

    % Linear interpolation to find exact crossing
    P1 = valid_powers(crossing_idx);
    P2 = valid_powers(crossing_idx + 1);
    alpha1 = valid_reals(crossing_idx);
    alpha2 = valid_reals(crossing_idx + 1);
    omega1 = valid_freqs(crossing_idx);
    omega2 = valid_freqs(crossing_idx + 1);

    % Interpolate power where alpha = 0
    if abs(alpha2 - alpha1) < 1e-12
        % Nearly parallel to axis, use midpoint
        P_boundary = (P1 + P2) / 2;
        omega_boundary = (omega1 + omega2) / 2;
        interp_quality = 0.5;
    else
        % Linear interpolation
        t = -alpha1 / (alpha2 - alpha1);  % Parameter where alpha = 0
        P_boundary = P1 + t * (P2 - P1);
        omega_boundary = omega1 + t * (omega2 - omega1);
        interp_quality = min(t, 1-t);  % Best when t ≈ 0.5
    end

    % Classify instability type
    if omega_boundary < 0.01
        instability_type = 'divergent';  % RBMK case: exponential growth
    else
        instability_type = 'hopf';  % Oscillatory instability
    end

    %% Build info structure
    info.eigenvalue_trace = eigenvalue_trace;
    info.frequency_trace = frequency_trace;
    info.crossing_index = crossing_idx;
    info.interpolation_quality = interp_quality;
    info.crossing_powers = [P1, P2];
    info.crossing_reals = [alpha1, alpha2];
    info.crossing_freqs = [omega1, omega2];
    info.detection_method = detection_method;
    info.zero_tolerance = zero_tolerance;
    info.instability_type = instability_type;

    %% Display results
    if opts.Verbose
        fprintf('STABILITY BOUNDARY DETECTION:\n');
        fprintf('  Detection method: %s\n', detection_method);
        fprintf('  Crossing found between %.0f MW and %.0f MW\n', P1, P2);
        fprintf('  Real parts: %+.6f → %+.6f\n', alpha1, alpha2);
        fprintf('  Frequencies: %.4f → %.4f rad/s\n\n', omega1, omega2);

        fprintf('STABILITY BOUNDARY POINT:\n');
        fprintf('  Power: %.1f MW\n', P_boundary);
        fprintf('  Instability type: %s\n', upper(instability_type));

        if strcmp(instability_type, 'hopf')
            fprintf('  Frequency: %.4f rad/s\n', omega_boundary);
            fprintf('  Period: %.2f seconds\n', 2*pi/omega_boundary);
            fprintf('  → Oscillatory instability (limit cycle)\n');
        else
            fprintf('  Frequency: %.4f rad/s (essentially zero)\n', omega_boundary);
            fprintf('  → Divergent instability (exponential growth)\n');
        end

        fprintf('  Interpolation quality: %.3f (1.0 = best)\n', interp_quality);

        if interp_quality < 0.2
            fprintf('  ⚠ Warning: Interpolation near edge, consider refining grid\n');
        end

        if strcmp(detection_method, 'near_zero') || strcmp(detection_method, 'closest_to_zero')
            fprintf('  ⚠ Note: Margins very small (|α| < %.0e), boundary may be approximate\n', zero_tolerance);
        end

        fprintf('================================\n\n');
    end
end

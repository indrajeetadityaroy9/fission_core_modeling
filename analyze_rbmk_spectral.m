function analyze_rbmk_spectral()
% ANALYZE_RBMK_SPECTRAL - High-Precision Stability Map
%
% UPDATED LOGIC:
%   - Force Dense Solver (eig) to ensure no roots are missed.
%   - Prioritize Growth Rate (Re > 0) over Frequency filtering.
%   - Detects slow spatial waves (f > 0.001 Hz).

    clear; clc;
    fprintf('============================================================\n');
    fprintf(' RBMK SPECTRAL STABILITY ANALYSIS (Fixed Tracking)\n');
    fprintf('============================================================\n');

    % 1. Load Parameters
    P = rbmk_parameters();
    p_base = P.accident; % Use the accident configuration

    % Ensure coupling allows spatial modes
    if p_base.Dn > 0.5
        p_base.Dn = 0.2;
    end

    fprintf('[CONFIG] Accident Mode (Flow=%.0f, Dn=%.1f)\n', p_base.m_flow, p_base.Dn);

    % 2. Simulation Setup
    powers_mw = linspace(200, 3200, 50);
    max_growth = zeros(size(powers_mw));
    osc_freq = zeros(size(powers_mw));

    fprintf('\nCalculating Spectra (N=12, Dense Solver)...\n');
    fprintf('%-10s | %-12s | %-12s | %-10s\n', 'Power(MW)', 'Growth Rate', 'Freq (Hz)', 'Status');
    fprintf('------------------------------------------------------------\n');

    for i = 1:length(powers_mw)
        target_mw = powers_mw(i);
        pow_frac = target_mw / 3200;

        % A. Equilibrium
        [y_eq, p_run] = compute_equilibrium(pow_frac, p_base, 'normal');

        % B. Jacobians
        [J0, Jtau] = get_jacobians(y_eq, p_run);

        % C. Spectral Eigenvalues (FORCE DENSE SOLVER)
        % Set large_system_threshold to Inf to force 'eig' instead of 'eigs'
        [eigs_spec, ~, ~] = chebyshev_eigenvalues(J0, Jtau, p_run.tau_flow, ...
            'N', 12, 'large_system_threshold', Inf);

        % D. INTELLIGENT SORTING
        % 1. Filter Trivial Null Modes (Translation Invariance)
        valid_eigs = eigs_spec(abs(eigs_spec) > 1e-4);

        if isempty(valid_eigs)
            dom_lambda = -1e-9;
        else
            % 2. Sort by Real Part (Descending) - Find the most unstable root
            [~, idx] = sort(real(valid_eigs), 'descend');
            sorted = valid_eigs(idx);
            dom_lambda = sorted(1);
        end

        % E. Store Results
        max_growth(i) = real(dom_lambda);

        % Only report frequency if it's a dynamic mode (not pure divergence)
        if abs(imag(dom_lambda)) > 1e-3
            osc_freq(i) = abs(imag(dom_lambda)) / (2*pi);
        else
            osc_freq(i) = 0; % Pure exponential growth/decay
        end

        if max_growth(i) > 1e-5
            status = 'UNSTABLE';
        else
            status = 'Stable';
        end

        % Log output
        if mod(i, 5) == 0 || i == 1 || i == length(powers_mw)
             fprintf('%6.0f     | %+8.4f     | %8.4f     | %s\n', ...
                target_mw, max_growth(i), osc_freq(i), status);
        end

        % Save Spectra
        if i == 1
            spec_low = eigs_spec; mw_low = target_mw;
        elseif i == length(powers_mw)
            spec_high = eigs_spec; mw_high = target_mw;
        end
    end

    %% 3. VISUALIZATION
    figure('Name', 'RBMK Stability Map', 'Color', 'w', 'Position', [100 100 1000 700]);

    % Plot 1: Bifurcation Diagram
    subplot(2, 2, [1 2]);
    plot(powers_mw, max_growth, 'b-', 'LineWidth', 2);
    hold on;
    yline(0, 'r--', 'Threshold');
    grid on;
    title('Stability Map: Growth Rate vs Power');
    ylabel('Growth Rate (Re(\lambda))');
    xlabel('Thermal Power (MW)');
    ylim([-0.5, 0.5]); % Zoomed out to see the real growth

    % Plot 2: Spectrum at Low Power
    subplot(2, 2, 3);
    plot_spectrum(spec_low, mw_low, 'r');

    % Plot 3: Spectrum at High Power
    subplot(2, 2, 4);
    plot_spectrum(spec_high, mw_high, 'b');
end

function plot_spectrum(eigs_data, mw, color)
    % Filter huge spurious roots
    view_mask = abs(real(eigs_data)) < 5 & abs(imag(eigs_data)) < 10;
    d = eigs_data(view_mask);

    plot(real(d), imag(d), 'o', 'Color', color, 'MarkerFaceColor', color);
    hold on; xline(0, 'k--'); grid on;

    % Highlight Unstable
    unstable = d(real(d) > 0);
    if ~isempty(unstable)
        plot(real(unstable), imag(unstable), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5);
    end

    title(sprintf('Spectrum @ %.0f MW', mw));
    xlabel('Growth'); ylabel('Freq');
    xlim([-1, 0.5]);
end

function [J0, Jtau] = get_jacobians(y_eq, p)
    % Central Difference
    n = length(y_eq);
    eps_pert = 1e-6 * max(abs(y_eq), 1e-6);
    J0=zeros(n); Jtau=zeros(n); Z_eq=y_eq;

    for j=1:n
        h=eps_pert(j);
        yp=y_eq; yp(j)=yp(j)+h; ym=y_eq; ym(j)=ym(j)-h;
        J0(:,j)=(rbmk_dynamics(0,yp,Z_eq,p)-rbmk_dynamics(0,ym,Z_eq,p))/(2*h);

        Zp=Z_eq; Zp(j)=Zp(j)+h; Zm=Z_eq; Zm(j)=Zm(j)-h;
        Jtau(:,j)=(rbmk_dynamics(0,y_eq,Zp,p)-rbmk_dynamics(0,y_eq,Zm,p))/(2*h);
    end
end

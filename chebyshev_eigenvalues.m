function [eigenvalues, augmented_matrix, info] = chebyshev_eigenvalues(J0, J_tau, tau, varargin)
% CHEBYSHEV_EIGENVALUES - Spectral Collocation DDE Eigenvalue Solver
%
% DESCRIPTION:
%   Computes the eigenvalues of a Delay Differential Equation (DDE) system
%   linearized around an equilibrium point using the Chebyshev spectral
%   collocation method.
%
% SYSTEM:
%   dy/dt = J0 * y(t) + J_tau * y(t-tau)
%
% INPUTS:
%   J0    - [n x n] Jacobian matrix with respect to current state y(t)
%   J_tau - [n x n] Jacobian matrix with respect to delayed state y(t-tau)
%   tau   - [scalar] Delay time
%
% PARAMETERS (Name-Value Pairs):
%   'N'          - [int] Number of Chebyshev intervals (default: 12).
%                  N=12 is usually sufficient for machine precision.
%   'eigs_count' - [int] Number of eigenvalues to return if using sparse solver.
%
% OUTPUTS:
%   eigenvalues      - [vector] Sorted complex eigenvalues (descending Real part).
%   augmented_matrix - [matrix] The constructed linear operator approximation.
%   info             - [struct] Diagnostic metadata.

    %% 1. Input Validation
    p = inputParser;
    addRequired(p, 'J0', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'J_tau', @(x) isnumeric(x) && ismatrix(x) && all(size(x)==size(J0)));
    addRequired(p, 'tau', @(x) isscalar(x) && x > 0);
    addParameter(p, 'N', 12, @(x) isscalar(x) && x >= 2);
    addParameter(p, 'eigs_count', 50, @(x) isscalar(x) && x > 0);
    addParameter(p, 'large_system_threshold', 200, @isscalar);

    parse(p, J0, J_tau, tau, varargin{:});

    N = p.Results.N;
    n = size(J0, 1);

    %% 2. Spectral Discretization
    % Get Chebyshev points and differentiation matrix on [-1, 1]
    [D, xi] = cheb_diff(N);

    % Map points from [-1, 1] to [-tau, 0]
    % theta = (tau/2) * (xi - 1)
    theta = (tau / 2) * (xi - 1);

    % Scale differentiation matrix D by d(xi)/d(theta) = 2/tau
    D_scaled = (2 / tau) * D;

    %% 3. Augmented Matrix Construction
    % The state vector U corresponds to [u_0; u_1; ...; u_N]
    % where u_k approx y(t + theta_k).
    % u_0 is current time (theta=0), u_N is delayed time (theta=-tau).

    dim_aug = n * (N + 1);
    A_aug = zeros(dim_aug, dim_aug);

    % Block 1: The DDE Boundary Condition (dy/dt = J0*y + Jtau*y_delayed)
    % Corresponds to the derivative at u_0
    A_aug(1:n, 1:n) = J0;                    % Coeff for u_0
    A_aug(1:n, end-n+1:end) = J_tau;         % Coeff for u_N

    % Blocks 2..N+1: Spectral Differentiation Constraints
    % du_k/dt = Sum(D_kj * u_j)
    for k = 1:N
        row_idx = (k * n) + 1 : (k + 1) * n;

        for j = 0:N
            col_idx = (j * n) + 1 : (j + 1) * n;
            A_aug(row_idx, col_idx) = D_scaled(k+1, j+1) * eye(n);
        end
    end

    %% 4. Solve Eigenvalues
    t_start = tic;

    if dim_aug < p.Results.large_system_threshold
        % Dense solver for smaller systems (Most robust)
        eigenvalues = eig(A_aug);
        solver_type = 'eig (dense)';
    else
        % Sparse solver for larger N
        k = min(p.Results.eigs_count, dim_aug - 2);
        try
            % Find eigenvalues with largest real part (most unstable)
            eigenvalues = eigs(sparse(A_aug), k, 'largestreal');
            solver_type = 'eigs (sparse)';
        catch
            warning('Sparse solver failed. Fallback to dense eig().');
            eigenvalues = eig(A_aug);
            solver_type = 'eig (dense fallback)';
        end
    end

    % Sort by Real part (descending)
    [~, sort_idx] = sort(real(eigenvalues), 'descend');
    eigenvalues = eigenvalues(sort_idx);

    %% 5. Output Packaging
    if nargout >= 2
        augmented_matrix = A_aug;
    else
        augmented_matrix = [];
    end

    info.N = N;
    info.tau = tau;
    info.dimension = dim_aug;
    info.computation_time = toc(t_start);
    info.solver = solver_type;
    info.chebyshev_points = theta;

    if ~isempty(eigenvalues)
        info.stability_margin = real(eigenvalues(1));
    else
        info.stability_margin = NaN;
    end
end

%% ========================================================================
%  Helper: Chebyshev Differentiation Matrix
%  Algorithm: Weideman and Reddy (2000)
%  ========================================================================
function [D, x] = cheb_diff(N)
    % 1. Chebyshev-Gauss-Lobatto points
    k = (0:N)';
    x = cos(pi * k / N);

    % 2. Barycentric weights
    c = [2; ones(N-1, 1); 2] .* (-1).^k;

    % 3. Construct D
    X = repmat(x, 1, N+1);
    dX = X - X';

    % Off-diagonal entries
    D = (c * (1 ./ c)') ./ (dX + eye(N+1));

    % Diagonal entries (row sum = 0)
    D = D - diag(sum(D, 2));
end

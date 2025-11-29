function sol = dde15s_new(ddefun,lags,history,tspan,options,varargin)
% DDE15S_NEW - Solves stiff Delay Differential Equations (DDEs)
%
% This is a modified version of dde15s that wraps MATLAB's built-in
% 'ode15s' (Variable Order, Variable Step) integrator to handle delays.
%
% CRITICAL UPDATE: This version includes the "solextend" patch required
% to run on MATLAB R2018b and newer versions. Without this patch,
% the solver will crash when concatenating solution structures.
%
% USAGE:
%   sol = dde15s_new(ddefun, lags, history, tspan, options, p)
%
% INPUTS:
%   ddefun   - Function handle @(t,y,Z) returning dy/dt
%   lags     - Vector of delay times (e.g., [2.0] for tau_flow)
%   history  - Initial history vector (y0) or function handle
%   tspan    - Integration interval [t0 tf]
%   options  - ODE options from odeset (e.g., RelTol, AbsTol)
%   varargin - Extra parameters passed to ddefun (e.g., p structure)
%
% OUTPUT:
%   sol      - Solution structure with fields:
%              .x - Time points
%              .y - State values (columns correspond to .x)
%              .discont - Discontinuity points (from Method of Steps)
%
% EXAMPLE (RBMK Simulation):
%   p = rbmk_parameters();
%   [y0, p, ~] = compute_equilibrium(0.5, p);
%   lags = p.tau_flow;
%   tspan = [0, 100];
%   opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
%   sol = dde15s_new(@rbmk_dynamics, lags, y0, tspan, opts, p);
%   plot(sol.x, sol.y(1,:));  % Plot n_L over time
%
% ALGORITHM:
%   Uses the Method of Steps to propagate discontinuities through the
%   delay interval. The solution is computed in segments of length <= tau,
%   where each segment uses the previous segment's solution as history.
%
% Author: Modified for RBMK Forensics Model (R2018b+ compatible)

% Check inputs
if nargin < 4
  error('dde15s_new:NotEnoughInputs', ...
        'Must specify ddefun, lags, history, and tspan.');
elseif nargin == 4
  options = [];
end

t0 = tspan(1);
tfinal = tspan(end);
if tfinal <= t0
  error('dde15s_new:InvalidTspan', 'Must have tspan(1) < tspan(end).')
end

% Initialize y0 from history
if isnumeric(history)
  temp = history;
else
  temp = feval(history, t0, varargin{:});
end
y0 = temp(:);

% Check for custom initial conditions
maxlevel = 4;
initialy = ddeget(options, 'InitialY', []);
if ~isempty(initialy)
  y0 = initialy(:);
  maxlevel = 5;
end

t = t0;
y = y0;

%% DISCONTINUITY HANDLING (Method of Steps)
% DDEs propagate discontinuities at t0, t0+tau, t0+2*tau, etc.
% We must integrate through each discontinuity explicitly.

if isempty(lags)
  discont = tfinal;
  minlag = Inf;
else
  lags = lags(:)';
  minlag = min(lags);
  if minlag <= 0
    error('dde15s_new:InvalidLags', 'All lags must be positive.')
  end

  % Start with initial discontinuity
  vl = t0;
  maxlag = max(lags);

  % Check for user-specified jumps (discontinuities in parameters)
  jumps = ddeget(options, 'Jumps', []);
  if ~isempty(jumps)
    indices = find(((t0 - maxlag) <= jumps) & (jumps <= tfinal));
    if ~isempty(indices)
      jumps = jumps(indices);
      vl = sort([vl jumps(:)']);
      maxlevel = 5;
    end
  end

  % Propagate discontinuities through multiple delay intervals
  discont = vl;
  for level = 2:maxlevel
    vlp1 = vl(1) + lags;
    for i = 2:length(vl)
      vlp1 = [vlp1 (vl(i) + lags)];
    end
    indices = find(vlp1 <= tfinal);
    vl = vlp1(indices);
    if isempty(vl)
      break;
    end
    nvl = length(vl);
    if nvl > 1
      vl = sort(vl);
      % Remove duplicate discontinuities (within numerical tolerance)
      indices = find(abs(diff(vl)) <= 10*eps*abs(vl(1:nvl-1))) + 1;
      vl(indices) = [];
    end
    discont = [discont vl];
  end

  % Sort and remove duplicates from final discontinuity list
  if length(discont) > 1
    discont = sort(discont);
    indices = find(abs(diff(discont)) <= 10*eps*abs(discont(1:end-1))) + 1;
    discont(indices) = [];
  end
end

% Ensure tfinal is included
if abs(tfinal - discont(end)) <= 10*eps*abs(tfinal)
  discont(end) = tfinal;
else
  discont = [discont tfinal];
end

% Remove points at or before t0
indices = find(discont <= t0);
discont(indices) = [];
discont = [t0 discont];

%% INITIALIZE SOLUTION STRUCTURE
sol.x = t0;
sol.y = y0;

%% MAIN INTEGRATION LOOP (Method of Steps)
% Integrate through each discontinuity interval
for i = 2:length(discont)
  distance = discont(i) - discont(i-1);
  nsteps = ceil(distance / minlag);
  stepsize = distance / nsteps;
  e = discont(i-1);

  for j = 1:nsteps
    b = e;
    e = b + stepsize;
    if j == nsteps
      e = discont(i);  % Ensure we hit the exact endpoint
    end

    % Call ode15s for the current step
    % The wrapper ddefcn computes Z (delayed states) from sol history
    solex = ode15s(@ddefcn, [b, e], sol.y(:,end), options, ...
                   ddefun, lags, history, sol, varargin{:});

    % Concatenate solutions
    if (i == 2) && (j == 1)
      sol = solex;
    else
      sol = solextend(sol, solex);
    end
  end
end

% Store discontinuity information
sol.discont = discont;

end  % End of main function


%% =========================================================================
%  SUBFUNCTION: solextend
%  Concatenates two solution structures from ode15s
%  INCLUDES CRITICAL PATCH FOR MATLAB R2018b+
%% =========================================================================
function solout = solextend(sol, solex)
% SOLEXTEND - Concatenate solution structures
%
% This function handles the internal data structures of ode15s, which
% changed in MATLAB R2018b to use 3D arrays with varying dimensions.
%
% --- R2018b PATCH EXPLANATION ---
% In older MATLAB versions, sol.idata.dif3d was a 2D array.
% In R2018b+, it became a 3D array where the 3rd dimension can vary.
% Direct concatenation [sol.idata.dif3d, solex.idata.dif3d] fails
% because the dimensions don't match.
%
% The fix: Pre-allocate a combined array with the maximum required
% dimensions, then copy data into the appropriate regions.

solout.solver = 'ode15s';
solout.x = [sol.x, solex.x(2:end)];
solout.y = [sol.y, solex.y(:, 2:end)];

% --- BEGIN PATCH FOR MODERN MATLAB (R2018b+) ---

% Concatenate kvec (1D vector - straightforward)
solout.idata.kvec = cat(2, sol.idata.kvec, solex.idata.kvec(2:end));

% Handle dif3d (3D array with potentially different dimensions)
[s1, s2, s3] = size(sol.idata.dif3d);
[~, s2n, s3n] = size(solex.idata.dif3d);

% Use superiorfloat to ensure consistent data type
dataType = superiorfloat(sol.x, solex.x);

% Pre-allocate combined array with maximum necessary dimensions
solout.idata.dif3d = zeros(s1, max(s2, s2n), s3 + s3n + 1 - 2, dataType);

% Insert old solution data
solout.idata.dif3d(:, 1:s2, 1:s3) = sol.idata.dif3d;

% Insert new solution data (skip first time point to avoid duplication)
solout.idata.dif3d(:, 1:s2n, s3+1:end) = solex.idata.dif3d(:, :, 2:end);

% --- END PATCH ---

solout.extdata = [];

end  % End of solextend


%% =========================================================================
%  SUBFUNCTION: ddefcn
%  Wrapper function that computes delayed states Z for the DDE
%% =========================================================================
function v = ddefcn(t, y, ddefun, lags, history, sol, varargin)
% DDEFCN - Compute dy/dt with delayed states
%
% This function is called by ode15s at each internal step.
% It evaluates the delayed states Z by:
%   1. Looking up y(t-tau) from the existing solution (sol)
%   2. Using the initial history for t-tau < t0
%
% INPUTS:
%   t        - Current time
%   y        - Current state vector
%   ddefun   - User's DDE function handle
%   lags     - Vector of delay times
%   history  - Initial conditions (vector or function)
%   sol      - Solution structure built so far
%   varargin - Extra parameters for ddefun

if isempty(lags)
  Z = [];
else
  nlags = length(lags);
  Z = zeros(length(y), nlags);

  for i = 1:nlags
    % Compute the delayed time
    tmlag = t - lags(i);

    % Clamp to available solution range (safety for numerical issues)
    tmlag = min(tmlag, sol.x(end));

    if tmlag <= sol.x(1)
      % Before solution start: use initial history
      if isnumeric(history)
        Z(:, i) = history(:);
      else
        Z(:, i) = feval(history, tmlag, varargin{:});
      end
    else
      % Within solution range: interpolate using deval
      Z(:, i) = deval(sol, tmlag);
    end
  end
end

% Call the user's DDE function
v = feval(ddefun, t, y, Z, varargin{:});

end  % End of ddefcn

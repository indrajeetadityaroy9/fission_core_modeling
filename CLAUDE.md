# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an RBMK nuclear reactor Hopf bifurcation analysis codebase that mathematically demonstrates why low-power operation was dangerous (the Chernobyl disaster mechanism). It models thermal-hydraulic oscillations arising from positive void coefficient + transport delay using a 16-state delay differential equation (DDE) system.

## Running the Analysis

**Main entry point:**
```matlab
main
```

This runs the Hopf bifurcation study (~5 minutes):
- Sweeps power from 100-2000 MW
- Computes steady states and eigenvalues via Chebyshev spectral method
- Identifies critical power threshold for oscillatory instability

**Output:** `rbmk_bifurcation_results.mat` containing all analysis data.

**Quick test at specific power level:**
```matlab
p = rbmk_parameters();
[y_ss, info] = compute_equilibrium(0.5, p);  % 50% power
[eigs, margin, info] = compute_stability(y_ss, p);
```

## Architecture

### Core Components

**Parameter Configuration:**
- `rbmk_parameters.m` - Returns 50+ calibrated parameters for the RBMK-1000 reactor model. All parameters are forensically calibrated from INSAG-7, OECD-NEA, and research literature. This is the single source of truth for reactor physics.

**Physics Model (16-state DDE):**
- `rbmk_dynamics.m` - Computes time derivatives for the two-region reactor model. State vector: `[n, C, α, T_f, T_m, I, X, c]` for lower and upper core regions (8 states each). Implements positive void coefficient, Doppler feedback, xenon dynamics, and SCRAM effects.

**Steady-State Solver:**
- `compute_equilibrium.m` - Finds equilibrium states via nonlinear least squares (lsqnonlin). Weighted residuals ensure xenon reaches true equilibrium.
- `compute_equilibrium_branch.m` - Computes steady states across power range with continuation method.

**Stability Analysis:**
- `compute_stability.m` - Linear stability analysis via Chebyshev pseudospectral method. Computes DDE eigenvalues by discretizing the delay interval with Chebyshev collocation.
- `chebyshev_eigenvalues.m` - Core spectral method implementation. Builds augmented system matrix and solves eigenvalue problem.
- `analyze_branch_stability.m` - Loops over steady states to compute eigenvalues at each power level.
- `find_hopf_point.m` - Detects where eigenvalues cross imaginary axis (Hopf bifurcation point).

**Bifurcation Sweep:**
- `hopf_analysis.m` - 2-step pipeline: compute steady-state branch → analyze stability via eigenvalues

**Visualization:**
- `visualize_hopf_spiral.m` - 3D spiral bifurcation diagram showing stability degradation
- `visualize_core_grid.m` - Core cross-section visualization

### Data Flow

```
rbmk_parameters
        ↓
compute_equilibrium_branch
        ↓
compute_equilibrium (per power level)
        ↓
analyze_branch_stability
        ↓
compute_stability → chebyshev_eigenvalues
        ↓
find_hopf_point → Hopf bifurcation point
```

## Key Physics

The model captures three RBMK design flaws:

1. **Positive void coefficient** (κ_V = +0.025): More steam → more reactivity → more power
2. **Graphite-tipped control rods** (ρ_tip = +0.005 = +1β): SCRAM initially adds positive reactivity
3. **Slow SCRAM insertion** (τ_c = 18s): Motor-driven rods take 18 seconds for full insertion

The 2-second coolant transport delay (`p.tau_flow`) creates a phase lag that enables Hopf bifurcation at low power (~762 MW).

## State Vector Convention

The 16-element state vector uses this ordering for each region (lower = indices 1-8, upper = 9-16):
1. `n` - Neutron density (normalized)
2. `C` - Delayed neutron precursors
3. `α` - Void fraction (steam)
4. `T_f` - Fuel temperature
5. `T_m` - Moderator temperature
6. `I` - Iodine-135 concentration
7. `X` - Xenon-135 concentration
8. `c` - Control rod position

Power is computed as `P(MW) = (n_L + n_U) * p.k_P` where `n_L = Y(1)` and `n_U = Y(9)` are normalized neutron densities.

## Common Modifications

**Adjust bifurcation sweep range:**
```matlab
[powers, eigs, margins] = hopf_analysis(p, ...
    'power_range', [100 2000], 'n_points', 50);
```

**Modify reactor parameters for sensitivity studies:**
```matlab
p = rbmk_parameters();
p.kappa_V = 0.020;  % Reduce void coefficient
```

**Run silent with custom tolerances:**
```matlab
[y_ss, info] = compute_equilibrium(0.5, p, 'Verbose', false, 'FunctionTol', 1e-10);
```

## Numerical Methods

- **Steady-state:** Levenberg-Marquardt via lsqnonlin with adaptive weighting for xenon equations (1000×) to ensure true equilibrium
- **Eigenvalues:** Chebyshev spectral collocation (N=12 points gives machine precision), augmented matrix dimension 208×208
- **Hopf detection:** Linear interpolation at stability margin zero-crossing

## Generated Output Files

After running `main`, the following files are created:
- `rbmk_bifurcation_results.mat` - All analysis data (power levels, eigenvalues, stability margins, Hopf point)
- `rbmk_analysis_*.log` - Detailed console log of the analysis session
- `hopf_sweep_details_*.txt` - Tab-separated eigenvalue data at each power level
- `steady_state_branch_*.txt` - Steady-state convergence data
- `stability_branch_*.txt` - Stability margin data

## Visualization (Post-Analysis)

After running `main`, generate visualizations from saved results:
```matlab
results = load('rbmk_bifurcation_results.mat');
visualize_hopf_spiral(results);  % 3D bifurcation diagram
```

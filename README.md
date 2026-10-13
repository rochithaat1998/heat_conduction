# 2D Transient Heat Conduction – FTCS (Explicit) and ADI (Implicit) in MATLAB


**Author:** Rochitha A T (AM24S012)

MATLAB codes that solve the 2D unsteady heat conduction equation in a rectangular bar using:

- **FTCS** – Forward Time Central Space (explicit)
- **ADI** – Alternating Direction Implicit (implicit, tridiagonal systems solved with the Thomas algorithm / TDMA)

The codes track the transient temperature field until steady state, plot mid-plane temperature profiles at selected times, report computation/CPU time, and a separate script computes the heat transfer through selected boundaries.

---

## Problem Definition

A long rectangular bar of length `L = 0.3 m` and height `H = 0.4 m` is initially at `T₀ = 0 °C`. Its four surfaces are then held at constant temperatures.

| Boundary | Location | Temperature |
|---|---|---|
| T1 | bottom (y = 0) | 40 °C |
| T2 | right (x = L) | 0 °C |
| T3 | top (y = H) | 10 °C |
| T4 | left (x = 0) | 0 °C |

| Property | Value |
|---|---|
| Thermal conductivity, k | 380 W/m·°C |
| Thermal diffusivity, α | 11.234 × 10⁻⁵ m²/s |

Grid: **31 × 41** points in x and y (Δx = Δy = 0.01 m).

### Tasks

1. Solve the transient temperature distribution using FTCS (Δt = 0.5, 0.1, 0.01, 0.001 s) and ADI (Δt = 0.5, 0.1 s). Plot temperature profiles at the mid-planes
   - x = 0.15 m (along y)
   - y = 0.2 m (along x)

   at t = 0, 0.1, 0.5, 1, 5, 10, 20, 40 s until steady state.
2. Compute the heat transfer for the y = 0, y = H and x = 0 planes.
3. Compare total computational time and CPU time for the explicit and implicit schemes.

**Convergence criterion:** the change in temperature between time levels must fall below 0.01 °C (the code checks the maximum nodal change, after all storage times have been passed).

---

## Governing Equation

```
∂T/∂t = α ( ∂²T/∂x² + ∂²T/∂y² )
```

Heat flux (Fourier's law): `q = -k ∇T`, with `qx = -k ∂T/∂x` and `qy = -k ∂T/∂y`.

---

## Numerical Methods

### FTCS (explicit)

```
T(i,j)^(n+1) = T(i,j)^n + rx (T(i+1,j) - 2T(i,j) + T(i-1,j))^n
                        + ry (T(i,j+1) - 2T(i,j) + T(i,j-1))^n

rx = αΔt/Δx²,   ry = αΔt/Δy²
```

**Stability:** `rx + ry ≤ 0.5`. With Δx = Δy = 0.01 m this gives Δt ≤ ≈ 0.222 s, so Δt = 0.5 s is unstable and was not run. Δt = 0.1, 0.01 and 0.001 s are stable.

### ADI (implicit)

Each time step is split into two sub-steps:

1. Implicit in x, explicit in y (row-wise tridiagonal systems)
2. Implicit in y, explicit in x (column-wise tridiagonal systems)

Each tridiagonal system is solved with a custom `TDMA` function (Thomas algorithm). The method is unconditionally stable, but Δt should still be commensurate with the grid spacing to capture gradients.

---

## Repository Structure

```
heat_conduction/
├── ftcs_heat_conduction.m        # FTCS (explicit) solver + plots + timing
├── adi_heat_conduction.m         # ADI (implicit) solver with TDMA + plots + timing
├── heat_transfer.m               # Steady-state solve + heat transfer through y=0, y=H, x=0
├── results/                      # Saved figures
└── README.md
```

> Rename the files above to match the names in your repository.

---

## How to Run

**Requirements:** MATLAB (the report was published with R2024a). No extra toolboxes are needed.

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/<your-repo>.git
   cd <your-repo>/heat_conduction
   ```
2. Open the script you want in MATLAB.
3. Set the time step at the top of the file:
   ```matlab
   dt = 0.001;   % FTCS: 0.1, 0.01, 0.001 (0.5 is unstable) | ADI: 0.1
   ```
4. Run the script. It produces:
   - Steady-state temperature contour
   - Temperature profiles along y at x = 0.15 m for the stored times
   - Temperature profiles along x at y = 0.2 m for the stored times
   - Total computational time (`tic/toc`) and CPU time (`cputime`) printed in the command window

> **ADI note:** `TDMA` is defined as a local function. In a script file, local functions must appear at the **end** of the file (MATLAB R2016b+). Move the `TDMA` definition to the bottom of the script if MATLAB reports an error.

### Heat transfer

Run `heat_transfer.m`. It solves the steady-state field with a Gauss–Seidel-type iteration, computes `qx`, `qy` with `gradient`, and integrates along the boundaries:

```matlab
Q_y0_plane = sum(qy(1,:)   * dx);   % y = 0
Q_yH_plane = sum(qy(end,:) * dx);   % y = H
Q_x0_plane = sum(qx(:,1)   * dy);   % x = 0
```

---

## Results

### Computation time (from the report)

| Scheme | Δt (s) | Total computation time (s) | CPU time (s) |
|---|---|---|---|
| FTCS | 0.001 | 1.5463 | 1.9062 |
| FTCS | 0.01 | 1.2751 | 1.8125 |
| FTCS | 0.1 | 1.2511 | 1.3281 |
| FTCS | 0.5 | Stability criterion not met | – |
| ADI | 0.1 | 1.0757 | 2.1562 |
| ADI | 0.5 | Not run (not commensurate with grid spacing) | – |

Timings depend on the machine and on each run; the values printed in the appendix of the report (e.g. 2.9769 s / 3.3750 s for FTCS at Δt = 0.001 s, 7.1845 s / 10.2969 s for ADI at Δt = 0.1 s) differ from the summary table, so re-run and record your own numbers.

### Heat transfer (steady state)

| Plane | Heat transfer |
|---|---|
| y = 0 | 63356.75 W |
| y = H | −15094.33 W |
| x = 0 | −56277.42 W |

The values are per unit depth of the bar (the problem is 2D), so strictly the units are W/m. The sign follows the direction of the +x/+y axes.



## Known Limitations / Notes

- **Profile plots:** the profile along y is plotted against `flip(y_dom)` while row 1 of `T` is y = 0 (T1 = 40 °C). As a result the plots show 40 °C at y = 0.4 m and 10 °C at y = 0, which is the reverse of the stated boundary conditions. Plotting against `y_dom` fixes this.
- **ADI coefficients:** the code uses `γ = αΔt/Δx²` in each sub-step, whereas the formulation in the report uses `r = αΔt/(2Δx²)` for each half step. Check that the two are consistent.
- **FTCS at Δt = 0.5 s** is not run because it violates the stability condition `rx + ry ≤ 0.5`. No stability check is coded; add one if you want the script to warn automatically.
- **Heat-transfer script grid:** it uses `Nx = 30`, `Ny = 40` while the transient solvers use 31 × 41, so the grids differ slightly. Boundary values from `gradient` use one-sided differences, and the corner singularities make the flux values grid-sensitive.
- A grid-independence study (optional in the assignment) was not performed.

---



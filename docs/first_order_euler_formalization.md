# First-order Euler solver: pseudo-code → formal implementation plan

This document mirrors the comments from `FirstorderEuler.cpp` and splits the work into practical parts for collaboration.

## Part 1 — Lock down data model with TriangularMesh as the single mesh parser

1. **Define the minimum solver state** (already scaffolded in `FirstorderEuler`):
   - `U[i] = [rho, rhou, rhov, rhoE]`
   - `Residual[i]`
   - `Area[i]`, `Perimeter[i]`
2. **Add strongly-typed edge records**:
   - interior edge: `(elemL, faceL, elemR, faceR, normal, length)`
   - boundary edge: `(elem, localFace, groupId, bcType, normal, length)`
3. **Use `TriangularMesh` to parse `.gri` once** and derive element/face connectivity, normals, lengths, and periodic pairing directly from that mesh object.
4. **Use `TriangularMesh` element area directly** for per-element area (`Area[e]`) and keep validation optional via `validateMeshOnLoad` for debugging.

## Part 2 — Primitive/conservative conversion utilities

Implement pure utility functions (unit-testable):

- `primitiveFromConserved(U, gamma)`
- `conservedFromPrimitive(rho, u, v, p, gamma)`
- `soundSpeed(rho, p, gamma)`
- pressure positivity guard

These are required by all boundary conditions and Riemann fluxes.

## Part 3 — Uniform initialization from inlet stagnation data

Physical meaning: start from a quiescent, spatially uniform guess derived from inlet total conditions and a chosen initial Mach number; this gives the solver a thermodynamically consistent first state before flux updates.

Use your current math in one function (`initUniformState`):

- Inputs: `M, alpha, rho0, a0, gamma, R`
- Compute:
  - `p0 = rho0 * a0^2 / gamma`
  - `T0 = a0^2 / (gamma R)`
  - `T = T0 / (1 + (gamma-1)/2 * M^2)`
  - `p = p0 * (T/T0)^(gamma/(gamma-1))`
  - `u = M a0 cos(alpha)`, `v = M a0 sin(alpha)`
  - `rho = p/(R T)`
  - `rhoE = p/(gamma-1) + 0.5 rho (u^2+v^2)`

Then broadcast to all cells.

## Part 4 — Numerical flux interface (collaborator integration point)

Create an abstract interface already used by your codebase style:

- `F = num_flux(UL, UR, n_hat)`
- plug-in options: Roe / HLLE / Rusanov

Boundary conditions should **only** provide ghost/right state `UR`; the flux function is shared.

## Part 5 — Boundary ghost-state builders

Per your comments, implement separate files:

- `inflow_ghost(UL, Tt, pt, alpha, n)`
- `outflow_ghost(UL, pout, n)`
- `wall_ghost(UL, n)`
- `periodic_ghost(UL, U_partner)`

Rule: edge loop asks BC module for `UR`, then calls the same `num_flux(UL, UR, n)`.

Implementation now encapsulates boundary-state logic in one pair of files:

- `include/solver/BoundaryConditions.h`
- `src/solver/BoundaryConditions.cpp`

This module provides boundary types and formulas for:

1. steady inflow (total conditions + Riemann invariant),
2. unsteady inflow wake modulation,
3. subsonic outflow from `J+` and entropy,
4. wall slip reflection,
5. periodic partner-state passthrough.

In the current `FirstorderEuler` scaffold, this is done inside
`computeEdgeFluxesAndWaveSpeeds()`, which is called every iteration by
`advanceToConvergedOrFinalTime()`. That is where boundary conditions are
imposed in the time-stepping loop.

For this mesh family, the boundary notation is read from `.gri` face titles
through `TriangularMesh` (`Curve1`, ..., `Curve8`) and then interpreted by a
small BC mapper:

- `Curve1`, `Curve5` → blade walls.
- `Curve2 ↔ Curve4` and `Curve6 ↔ Curve8` are periodic pairs in `.gri`; these
  are converted to interior-like interfaces by mesh preprocessing (not treated
  by the boundary ghost-state function).
- Remaining non-periodic curves are mapped to inflow/outflow policy.
- Curve3/Curve7 roles are configurable in `SolverConfig` (`inflowCurve`, `outflowCurve`).

## Part 6 — Single edge loop for residual + wave-speed

Implement one kernel:

`compute_residual_and_wavespeed(mesh, U, bc, flux) -> (R, S)`

For each edge:

1. build `(UL, UR)`
2. compute edge flux `Fhat`
3. accumulate residual to owner/neighbor cells
4. compute spectral radius `lambda = |u.n| + c`
5. accumulate `S[cell] += lambda * edge_length`

This keeps CFL logic and residual assembly consistent.

## Part 7 — Time stepping (global vs local)

- Global: `dt = CFL * min(d_i) / avg(s_i)`
- Local: `dt_i = 2 * Area_i * CFL / S_i`
- Update: `U_i <- U_i - dt_i / Area_i * R_i`

Track `||R||_2` every N iterations and stop at tolerance.

## Part 8 — Post-processing hooks

After convergence:

- cell Mach number
- entropy
- blade `cp`
- force coefficients `(cx, cy)` using your normalization

Keep post-processing decoupled from the solver kernel for maintainability.

## Suggested teammate split

- **Teammate A**: mesh/connectivity readers + validation
- **Teammate B**: BC ghost-state functions
- **Teammate C**: numerical flux implementations
- **You (integration)**: edge-kernel + time-marching driver + convergence output



## Why use `TriangularMesh` as the single geometry source?

Short answer: yes, this is simpler and keeps the Euler solver focused.

- `TriangularMesh` already gives nodes, elements, faces, normals, lengths, and periodic relationships from `.gri`.
- `FirstorderEuler` now builds `interiorFaces_` and `boundaryFaces_` from that object instead of manually parsing `I2E/B2E/In/Bn/periodicEdges`.
- We now use `TriangularMesh::area(i)` directly, same as faces/normals/connectivity, to keep mesh data flow uniform.
- The marching loop still uses compact in-memory arrays for speed and clarity.

So the workflow becomes: **read `.gri` via `TriangularMesh` once + build solver arrays once**.


## Run modes now exposed in `FirstorderEuler`

To support team split with first-order vs second-order work, the solver now exposes explicit run-mode entry points:

- `runSteadyGlobal()`
- `runSteadyLocal()`
- `runUnsteadyGlobal()`

The default `main.cpp` accepts:

- `--mode steady-global`
- `--mode steady-local`
- `--mode unsteady-global`

with `--mesh <file.gri>` (or `--prefix <mesh_prefix>`, which maps to `<mesh_prefix>.gri` in `main.cpp`).

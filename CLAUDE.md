# LilyPad.jl

## Agent Routing Note
- This repo's local CLAUDE file is this project-root `CLAUDE.md`.
- Read the global `~/.claude/CLAUDE.md` first, then this local file.
- Write top-level persistent directives to the global `~/.claude/CLAUDE.md` file.
- On conversation compaction or resume, reread both the global file and this project-root `CLAUDE.md` before continuing.


## Package Purpose
This package extends the WaterLily.jl package with a semi-implicit time stepping scheme to increase simulation speed for interactive simulations. 

## Key Types

- `LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}` — thin wrapper around `WaterLily.Flow`; overrides `mom_step!`, `CFL`, and `BDIM!` with a semi-Lagrangian BDF2 step; holds two velocity buffers `u⁺`, `ū`; all other field access forwards to the inner `flow` via `getproperty`/`setproperty!`
- `LilyPadSim` — plain function (constructor alias), not a struct; returns `WaterLily.Simulation` built with `flow_ctor=LilyFlow`
- `WaterLily.AbstractPoisson` / `MultiLevelPoisson` — reused unchanged for pressure projection

## Entry Points

- `LilyPadSim(dims, uBC, L; Δt=1.5, kwargs...)` — same kwargs as `WaterLily.Simulation`; `Δt` sets the fixed timestep
- `sim_step!(sim[, t_end]; kwargs...)` — advance simulation (dispatches through `WaterLily.Simulation`)

## Source Layout

- `src/LilyPad.jl` — module entry point; `LilyPadSim` constructor alias, includes, exports
- `src/LilyFlow.jl` — `LilyFlow` struct + `getproperty`/`setproperty!` forwarding + `mom_step!` / `CFL` / `BDIM!` overrides
- `src/SemiAdvect.jl` — `departure` (RK2) and the two-level `advect!`
- `src/util.jl` — `quad` interpolation helpers
- `test/runtests.jl` — backend dispatch (mirrors WaterLily)
- `test/maintests.jl` — @testsets for each source file

## Environment Notes
A legacy code is in the ../lily-pad/ folder. The WaterLily source is in ../WaterLily(.jl)/.
WaterLily version: 1.10 (uuid: ed894a53-35f9-47f1-b17f-85db9237eebd); `vinterp` needs its component-indexed `_interp(x, u, i)` (WaterLily.jl#371)

## Critical Comparison Guardrail (Do Not Skip)
- **Always match viscosity between WL and LP before any timing or solver comparison.**
- **For inviscid comparisons, set `ν=0` explicitly for both solvers. Do not rely on defaults.**
- If using `Re`-based setup, confirm both cases map to the same `ν` numerically before running.
- Every benchmark/profiling script should print `ν` in its header and report it with results.
- Any result produced with mismatched `ν` is invalid and must be rerun.

## Development Order (Required)
Always execute work in this order:
1. Algorithm validation first
2. Implementation decisions second
3. Correctness testing third
4. Benchmarking last

Benchmarking is blocked until the first three stages are explicitly checked off.
Do not use performance results to justify an algorithm or method variant that has
not yet been considered and discussed.

## Testing Status
**Passing tests do NOT validate algorithm correctness.** Tests exist for:
- SemiAdvect: departure RK2 structure, scalar component advection on staggered faces, boundary-adjacent transport checks
- LilyFlow: one-step uniform-flow fixed-point behavior, divergence finiteness after projection, second-step regression (guards against scale_u! bug)
- LilyPadSim: one-step wrapper parity with WaterLily `Simulation`, timestep progression, second-step regression, periodic Taylor-Green steadiness
- LilyPadBiotExt: one-step parity with `BiotSimulation`, `symmetry` and `perdir` pass-through

Test suite runs on `Array` and/or `CuArray`, selected like WaterLily's `WATERLILY_BACKENDS`: `LILYPAD_BACKENDS=cpu|cuda|all` (default `cpu`). CUDA is not a test dependency; it is installed on demand when requested and skipped with a warning if not functional. Run GPU tests locally with `LILYPAD_BACKENDS=all` before any PR touching kernels or tests.
- GPU-safe test closures: capture only isbits values (use `f0` literals, not `T(...)`), and no `x[1:2]`-style slices of `SVector`s (they allocate).

## Scheme Notes (Semi-Lagrangian BDF2)
- `(3uⁿ⁺¹ - 4uⁿ(x) + uⁿ⁻¹(x⁻))/2Δt = -∇pⁿ⁺¹ + ν∇²u`, one projection per step (Xiu et al. 2005); `x`, `x⁻` backtrack one step each under `u⁺=(3uⁿ-uⁿ⁻¹)/2` and `ū=(uⁿ+uⁿ⁻¹)/2`
- `mom_project!` takes `w::Int`, so `p` is scaled by 2/3 around a `w=1` projection
- No CFL constraint from advection; `CFL(flow::LilyFlow)` returns `Δt[end]` unchanged to keep the step fixed
- `Δt` set at construction time is the user knob for speed/accuracy tradeoff; larger `Δt` is faster but less accurate
- Domain BCs are applied in WL's `mom_project!`, not in the advection kernel
- Viscosity diffuses each BDF2 sample over its age (`quad` variance `2νΔt` and `4νΔt`); stable for `νΔt ≤ 1/3`; `g` and `udf` throw
- Periodic `perdir`: interpolation wraps instead of clamping, and the face at index 2 is advected
- On WL bumps, diff the overridden functions against upstream, not just their signatures

## Timestep Control
- `Δt` is passed through `LilyPadSim(...; Δt=...)` → `WaterLily.Simulation` → `LilyFlow` → `Flow`
- `CFL(a::LilyFlow)` returns `a.Δt[end]` unchanged so `mom_step!`'s single `push!(Δt, CFL(a))` keeps the step fixed
- There is no adaptive mode — the semi-Lagrangian scheme always runs at the user-supplied `Δt`

## Legacy Timestep Reference
- Legacy repo: `../lily-pad/LilyPad/BDIM.pde`
- `dt == 0` means adaptive mode via `setDt()` and `checkCFL()`
- Main legacy circle demo: `../lily-pad/LilyPad/LilyPad.pde`
- The legacy circle demo constructs `BDIM(n,n,1.5,body)`, so the reference semi-Lagrangian circle run uses a constant `dt = 1.5`

## Current Status & Blockers
- **README case** (`examples/readme_circle_compare.jl`: circle L=32, Re=500, ν=0.064 for both, LP `Δt=2.5`, peak-to-peak over t=30–60, RTX A2000 Laptop GPU):

	|           | Cd    | Cd′   | CL    | St    | s/(tU/L) |
	|-----------|-------|-------|-------|-------|----------|
	| WaterLily | 1.742 | 0.128 | 1.410 | 0.256 | 0.271    |
	| LilyPad   | 1.913 | 0.148 | 1.564 | 0.254 | 0.020    |

	- LP vs WL: Cd +10%, Cd′ +16%, CL amplitude +11%, St −1%.
	- Speed-up 11–13×: 7× from `Δt` (WL's viscous CFL gives `Δt≈0.35`), 1.5× per step since the 2D case is kernel-launch bound (GPU busy 20–27%) and LP launches 161 vs 248 kernels per step.
	- Always average over whole cycles. Lift mean is not reported: its cycle-to-cycle standard error is too large to resolve.
- **3D sphere, GPU, ν=0, LP `Δt=2.5`** (s/(tU/L)): L=32 WL 2.60, LP 0.77 (3.4×); L=64 WL 50.3, LP 12.5 (4.0×). Large 3D grids are compute bound: LP costs 1.26× WL per step; advection ~50% and projection ~38% of LP GPU time; 3 MG iterations per projection vs WL's 1.
- **Viscous validation**: Taylor–Green decay error ≤0.2% (same as WL) for `νΔt ≤ 0.15`, effective ν +1.9% at `νΔt = 1/3`. Circle L=16, Re=200 vs WL: Cd +3.5%, CL amplitude +4%, St +1%.
- **Inviscid (ν=0) wakes collapse to grid scale** at every L; compare schemes only with ν>0.
- **Allocation note**: the README circle case still shows nontrivial allocations during body remeasurement in both WL and LP. This is likely a separate type-stability or body-measurement issue and is deferred for later investigation.
- **Current allocation-testing rule**: measure after initialization with `remeasure=false` and compare allocations per unit simulated time (`tU/L`), not per step.
- **Constraint**: treat performance numbers as provisional until each stage's correctness checks pass.

## Active Work
1. ✅ Refactored to WL `Flow-refactor` branch: `LilyFlow <: AbstractFlow{D,T}` thin wrapper, `LilyPadSim` constructor alias, `SemiMomStep.jl` deleted, `fixed_dt` removed.
2. ✅ GPU tests enabled: CUDA added as optional test dependency; closures in tests use `f0` literals to be GPU-compatible; `Pkg.test()` runs on both `Array` and `CuArray`.
3. ✅ Ensure working interoperability with BiotSavartBCs.jl capability for "unbounded" simulations.
4. ✅ Viscous diffusion along the trajectories (`νΔt ≤ 1/3`); WL vs LP benchmarks updated.
6. Turbulence: eddy viscosity in the trajectory diffusion, Lagrangian dynamic model, monotone `quad` limiter.
7. Trajectory-crossing time-step limit from `WaterLily.S`, at finite Re.
5. ✅ Develop real-time interactive 2D fluid simulation and vizualization example using ParticleViz and PixelBodies.
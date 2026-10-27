# LilyPad.jl

## Agent Routing Note
- This repo's local CLAUDE file is this project-root `CLAUDE.md`.
- Read the global `~/.claude/CLAUDE.md` first, then this local file.
- Write top-level persistent directives to the global `~/.claude/CLAUDE.md` file.
- On conversation compaction or resume, reread both the global file and this project-root `CLAUDE.md` before continuing.


## Package Purpose
This package extends the WaterLily.jl package with a semi-implicit time stepping scheme to increase simulation speed for interactive simulations. 

## Key Types

- `LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}` — thin wrapper around `WaterLily.Flow`; overrides `mom_predict!`, `mom_correct!`, `CFL`, and `BDIM!` with semi-Lagrangian logic; all field access forwards to the inner `flow` via `getproperty`/`setproperty!`
- `LilyPadSim` — plain function (constructor alias), not a struct; returns `WaterLily.Simulation` built with `flow_ctor=LilyFlow`
- `WaterLily.AbstractPoisson` / `MultiLevelPoisson` — reused unchanged for pressure projection

## Entry Points

- `LilyPadSim(dims, uBC, L; Δt=1.5, kwargs...)` — same kwargs as `WaterLily.Simulation`; `Δt` sets the fixed timestep
- `sim_step!(sim[, t_end]; kwargs...)` — advance simulation (dispatches through `WaterLily.Simulation`)

## Source Layout

- `src/LilyPad.jl` — module entry point; `LilyPadSim` constructor alias, includes, exports
- `src/LilyFlow.jl` — `LilyFlow` struct + `getproperty`/`setproperty!` forwarding + `mom_predict!` / `mom_correct!` / `CFL` / `BDIM!` overrides
- `src/SemiAdvect.jl` — `departure` (RK2, linear `vinterp` per face) and `advect!` (quadratic `quad` sample at the departure point, with viscous smoothing)
- `src/util.jl` — `vinterp` and `quad` interpolation helpers
- `ext/LilyPadBiotExt.jl` — `LilyBiotSim` for BiotSavartBCs domains
- `examples/readme_circle_compare.jl` — README circle case (WL vs LP, forces, vorticity, timing); `examples/biot_circle_gpu.jl` — BSBC circle
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

## Scheme Notes (Semi-Lagrangian two-step)
- Predictor: transport `uⁿ` along `uⁿ` over `Δt`, BDIM, project. Corrector: transport `uⁿ - (Δt/2)μ₀∇p` along `(uⁿ+u_pred)/2`, BDIM, project. Two projections per step, like WL
- Viscosity diffuses each departure sample with `quad` variance `2νΔt` (predictor and corrector); weights stay positive for `νΔt ≤ 1/3`, which the constructor checks
- No CFL constraint from advection; `CFL(flow::LilyFlow)` returns `Δt[end]` unchanged to keep the step fixed
- `Δt` set at construction time is the user knob for speed/accuracy tradeoff; larger `Δt` is faster but less accurate
- Domain BCs are applied in WL's `mom_project!`, not in `mom_predict!`/`mom_correct!`
- `g` and `udf` throw
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
- **README case** (`examples/readme_circle_compare.jl gpu`: circle L=64, Re=500, ν=0.128 for all, whole lift cycles (upward zero crossings) in t=30–60, RTX A2000 Laptop GPU, best of 3 timing windows):

	|                | Cd    | Cd′   | CL    | St    | s/(tU/L) |
	|----------------|-------|-------|-------|-------|----------|
	| WaterLily      | 1.724 | 0.130 | 1.431 | 0.258 | 1.373    |
	| LilyPad Δt=1.5 | 1.663 | 0.131 | 1.469 | 0.259 | 0.344    |
	| LilyPad Δt=2.5 | 1.644 | 0.131 | 1.525 | 0.259 | 0.202    |

	- Cd drops and CL rises with Δt. The WL time is slower than the 0.799 s measured on the advection-reflection branch with identical forces, so the GPU was likely throttled; compare speed-ups within one run only.
	- `νΔt ≤ 1/3` caps this case at Δt≈2.6.
	- Always average over whole cycles. Lift mean is not reported: its cycle-to-cycle standard error is too large to resolve.
- **Advection-reflection** (PR #18, branch `adv-reflection`, worktree `../LilyPad-ar`) is open and unmerged pending further study. On the README case it is closer to WL (Δt=2.5: Cd 1.741, CL 1.422) and allows `νΔt ≤ 2/3`.
- **Open question**: both schemes carry about 2× WL's near-body grid-scale vorticity at high Re; source not found. Uncommitted notes, scripts and data are in `examples/sphere_dt_investigation.md` (start there), `examples/pareto/`, `examples/sphere_dt/`, hidden via `.git/info/exclude`.
- **Allocation note**: the README circle case still shows nontrivial allocations during body remeasurement in both WL and LP. This is likely a separate type-stability or body-measurement issue and is deferred for later investigation.
- **Current allocation-testing rule**: measure after initialization with `remeasure=false` and compare allocations per unit simulated time (`tU/L`), not per step.
- **Constraint**: treat performance numbers as provisional until each stage's correctness checks pass.

## Active Work
1. ✅ Refactored to WL `Flow-refactor` branch: `LilyFlow <: AbstractFlow{D,T}` thin wrapper, `LilyPadSim` constructor alias, `SemiMomStep.jl` deleted, `fixed_dt` removed.
2. ✅ GPU tests enabled: CUDA added as optional test dependency; closures in tests use `f0` literals to be GPU-compatible; `Pkg.test()` runs on both `Array` and `CuArray`.
3. ✅ Ensure working interoperability with BiotSavartBCs.jl capability for "unbounded" simulations.
4. ✅ Viscous diffusion along the trajectories (`νΔt ≤ 1/3`, PR #20); README benchmarks updated.
5. ✅ Develop real-time interactive 2D fluid simulation and vizualization example using ParticleViz and PixelBodies.
6. Advection-reflection (PR #18): study before merging; near-body grid-scale noise at high Re.

# LilyPad.jl

## Agent Routing Note
- This repo's local CLAUDE file is this project-root `CLAUDE.md`.
- Read the global `~/.claude/CLAUDE.md` first, then this local file.
- Write top-level persistent directives to the global `~/.claude/CLAUDE.md` file.
- On conversation compaction or resume, reread both the global file and this project-root `CLAUDE.md` before continuing.

## Communication
- Use plain language. Don't invent labels, abbreviations or jargon ("H1", "cut", "the split"); name the thing itself.
- Only use notation that already appears in the code, the notes, or the user's own messages, and define any symbol the first time it's used.
- Keep answers short and direct: what you think, why, and what you'd do next.

## Package Purpose
This package extends the WaterLily.jl package with a semi-implicit time stepping scheme to increase simulation speed for interactive simulations. 

## Key Types

- `LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}` — thin wrapper around `WaterLily.Flow`; overrides `mom_step!`, `CFL`, and `BDIM!` with a semi-Lagrangian advection-reflection step; holds the reflected field `r`, the last final-projection correction `δp` (initial guesses only) and the cell-centred backtrack `disp`; all other field access forwards to the inner `flow` via `getproperty`/`setproperty!`
- `LilyPadSim` — plain function (constructor alias), not a struct; returns `WaterLily.Simulation` built with `flow_ctor=LilyFlow`
- `WaterLily.AbstractPoisson` / `MultiLevelPoisson` — reused unchanged for pressure projection

## Entry Points

- `LilyPadSim(dims, uBC, L; Δt=1.5, kwargs...)` — same kwargs as `WaterLily.Simulation`; `Δt` sets the fixed timestep
- `sim_step!(sim[, t_end]; kwargs...)` — advance simulation (dispatches through `WaterLily.Simulation`)

## Source Layout

- `src/LilyPad.jl` — module entry point; `LilyPadSim` constructor alias, includes, exports
- `src/LilyFlow.jl` — `LilyFlow` struct + `getproperty`/`setproperty!` forwarding + `mom_step!` / `CFL` / `BDIM!` overrides
- `src/SemiAdvect.jl` — `departure` (RK2, quadratic at the midpoint) and `advect!` (one backtrack per cell centre, cubic to the faces)
- `src/util.jl` — `vinterp`, `vquad` and `quad` interpolation helpers
- `examples/readme_circle_compare.jl`, `examples/tgv_dissipation.jl` — README circle and TGV Re=1600 harness
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
- LilyPadBiotExt: `symmetry` and `perdir` pass-through, Biot-Savart BCs slow the inflow face (< 0.97U), every projection converges (V-cycles < `itmx`). No accuracy comparison with `BiotSimulation`: after an impulsive start the old one-step parity only matched WL's first-step error
- Dependencies: `[sources]` points WaterLily at PR #245 (`poisson_kwargs`, `mom_project!` keywords) and BiotSavartBCs at its `poisson-kwargs` branch until both are released; after changing a git source run `Pkg.update("WaterLily")`, since `Pkg.resolve` keeps the old tree

Test suite runs on `Array` and/or `CuArray`, selected like WaterLily's `WATERLILY_BACKENDS`: `LILYPAD_BACKENDS=cpu|cuda|all` (default `cpu`). CUDA is not a test dependency; it is installed on demand when requested and skipped with a warning if not functional. Run GPU tests locally with `LILYPAD_BACKENDS=all` before any PR touching kernels or tests.
- GPU-safe test closures: capture only isbits values (use `f0` literals, not `T(...)`), and no `x[1:2]`-style slices of `SVector`s (they allocate).

## Scheme Notes (Semi-Lagrangian advection-reflection)
- Step (Zehnder, Narain & Thomaszewski 2018): advect `uⁿ` over `Δt/2` under `(5uⁿ-uⁿ⁻¹)/4`, BDIM, `mom_project!(w=2)`; reflect `r = u - (Δt/2)μ₀∇p₁` (= `2Pf - f`); advect `r` over `Δt/2` under `(3u^{n+1/2}-uⁿ)/2`, BDIM, add `Δtμ₀∇p₁` back, `mom_project!(w=1)`, which then solves for the physical pressure. Two standard projections per step for any `AbstractPoisson` (same count as WL's predictor-corrector); `δp` makes both initial guesses close (2.06 V-cycles/step at Δt=2.5, 3.11 at Δt=5)
- Why: a projection removes the gradient part of the transported field; at large Δt during vortex-sheet roll-up that removal was the excess dissipation (multistep schemes extrapolate the parcel's pressure history, which fails when it changes within a step). Reflection is an isometry, so only the final projection, whose input is nearly solenoidal, removes energy. The pressure impulse acts at the trajectory midpoint with a fresh projection every step, so there is no explicit-pressure feedback
- Backtracks are computed once per cell centre (RK2, quadratic interpolant at the midpoint) and interpolated to the faces with a 4-point cubic (linear where the stencil leaves the interior). Linear interpolation of the backtrack velocity biases it by ~h²∂²u/8 once the Courant number exceeds ~1 (1–3% of U in 2–3-cell shear layers); the quadratic midpoint removes that. Per cell: one backtrack instead of D; 63 vs 88 ms/step at 128³. `quad` loops over its 3^D stencil: unrolling it was 40% slower in 3D
- `reflect` evaluates the correction on ghost faces at the cell `BC!` would copy from (zero-gradient, periodic, zero on normal faces), so `r` keeps `u`'s domain BCs without another `BC!`; stale ghosts cost 30% accuracy on slip walls
- Multistep (BDF2/BDF3/blend, BDF4) variants were tried first; full BDF3 runs oscillate and gain energy after roll-up
- Explicit pressure sampled off the arrival point from the previous step (trajectory midpoint, SL Crank-Nicolson, path corrections) is unstable: pressure errors feed back with gain `|1-g|`
- Covector transport `P[Jᵀu(x_d)]` blows up in 3D TGV even at Δt=0.25 with `J` from neighbouring departure points, so it is not the RK4/ForwardDiff map
- ν=0 at small Δt gains energy (64³ TGV, Δt=1 E(20)=0.053; Δt=0.5 blows up): the semi-Lagrangian transport's own issue (old blend at Δt=0.5 too); reflection halves the per-transport interpolation damping
- No CFL constraint from advection; `CFL(flow::LilyFlow)` returns `Δt[end]` unchanged to keep the step fixed
- `Δt` set at construction time is the user knob for speed/accuracy tradeoff; larger `Δt` is faster but less accurate
- Domain BCs are applied in WL's `mom_project!`, not in the advection kernel
- Viscosity diffuses each half-step sample with `quad` variance `2ν(Δt/2)`; weights stay positive for `νΔt ≤ 2/3`, which the constructor checks; `g` and `udf` throw
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
- **README case** (`examples/readme_circle_compare.jl`: circle L=64, Re=500, ν=0.128 for all, whole lift cycles (upward zero crossings) in t=30–60, RTX A2000 Laptop GPU, best of 3 timing windows):

	|                | Cd    | Cd′   | CL    | St    | s/(tU/L) |
	|----------------|-------|-------|-------|-------|----------|
	| WaterLily      | 1.724 | 0.130 | 1.431 | 0.258 | 0.799    |
	| LilyPad Δt=2.5 | 1.741 | 0.129 | 1.422 | 0.260 | 0.090    |
	| LilyPad Δt=5   | 1.766 | 0.129 | 1.438 | 0.260 | 0.051    |

	- L=32 (old README case): WL Cd 1.747, CL 1.421, St 0.256; LP Δt=2.5 1.716, 1.389, 0.263 (St gap is LP resolution error, gone at L=64); old blend Δt=2.5 1.847, 1.497, 0.257.
	- Two standard projections per step; at 128³ a step costs ~63 ms at Δt=5 (old blend ~53 ms).
	- Always average over whole cycles. Lift mean is not reported: its cycle-to-cycle standard error is too large to resolve.
- **Resolved circle D=64, Re=200** (ν=0.32 for both, zero-crossing St): WL Cd 1.700, Cd′ 0.057, CL 0.990, St 0.2417; LP Δt=2 Cd 1.740, Cd′ 0.057, CL 1.003, St 0.2418 (9× faster). D=64, Re=500 (ν=0.128): WL Cd 1.724, Cd′ 0.130, CL 1.431, St 0.258; LP Δt=2.5 1.742, 0.129, 1.422, 0.260; Δt=5 1.766, 0.129, 1.436, 0.260 (8.9×, 17.8× faster). The L=32 St gap (+2.7%) shrinks to +0.8% at D=64 for both Δt, so it is LilyPad's resolution error at L=32, not the time step. Peak counting double-counts WL's small-step lift wiggles (St 0.447); use zero crossings.
- **TGV Re=1600, 128³** (DNS peak `-dE/dt` 0.0129 at t=9.0, E(20)=0.0215), rms `-dE/dt` error vs DNS on t∈[2,12]: reflection Δt=1 0.00026, Δt=2.5 0.00015, Δt=5 0.00026, Δt=7.5 0.00073, Δt=10 0.00140; old blend 0.00021, 0.00038, 0.00138, 0.00249, 0.00272. The Δt=5 roll-up bump is gone. Cost to t=20 (`examples/tgv_dissipation.jl`): WL 2.05–2.48 s/(tU/L) (run to run), LP Δt=1 2.54, 2.5 0.98, 5 0.50, 7.5 0.34; peak ν⟨ω²⟩ WL 0.0089, LP Δt=1 0.0115, 2.5 0.0107, 5 0.0100, 7.5 0.0087 (DNS 0.0128).
- **3D sphere L=32 (192×128×128), Re=500, ν=0.064 for both, GPU**: Cd WL 0.524, LP Δt=2.5 0.527, Δt=5 0.548; s/(tU/L) WL 3.285, LP 1.043 (3.2×), 0.525 (6.3×). (Old BDF2 ν=0 numbers: L=32 3.4×, L=64 4.0×; large 3D grids are compute bound, a reflection step costs about the same as a blend step.)
- **Viscous validation**: Taylor–Green decay error ≤0.2% (same as WL) for `νΔt ≤ 0.15`, effective ν +1.4% at `νΔt = 1/3` (blend). Circle L=16, Re=200 vs WL (BDF2 scheme): Cd +3.5%, CL amplitude +4%, St +1%.
- **Inviscid (ν=0) wakes collapse to grid scale** at every L; compare schemes only with ν>0.
- **Allocation note**: the README circle case still shows nontrivial allocations during body remeasurement in both WL and LP. This is likely a separate type-stability or body-measurement issue and is deferred for later investigation.
- **Current allocation-testing rule**: measure after initialization with `remeasure=false` and compare allocations per unit simulated time (`tU/L`), not per step.
- **Constraint**: treat performance numbers as provisional until each stage's correctness checks pass.

## Active Work
1. ✅ Refactored to WL `Flow-refactor` branch: `LilyFlow <: AbstractFlow{D,T}` thin wrapper, `LilyPadSim` constructor alias, `SemiMomStep.jl` deleted, `fixed_dt` removed.
2. ✅ GPU tests enabled: CUDA added as optional test dependency; closures in tests use `f0` literals to be GPU-compatible; `Pkg.test()` runs on both `Array` and `CuArray`.
3. ✅ Ensure working interoperability with BiotSavartBCs.jl capability for "unbounded" simulations.
4. ✅ Viscous diffusion along the trajectories (`νΔt ≤ 1/3`); WL vs LP benchmarks updated.
5. ✅ Develop real-time interactive 2D fluid simulation and vizualization example using ParticleViz and PixelBodies.
6. ✅ Turbulence: Δt≥5 roll-up over-dissipation fixed by advection-reflection with cell-centred quadratic backtracks Open: ν=0 small-Δt energy gain.
7. Trajectory-crossing time-step limit from `WaterLily.S`, at finite Re.
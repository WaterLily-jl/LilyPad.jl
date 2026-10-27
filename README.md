# LilyPad.jl

LilyPad.jl is a semi-Lagrangian momentum-step variant of WaterLily.jl.

The API is intentionally the same at the simulation level, so users can switch
the method by changing the constructor only.

## Method Difference

- WaterLily: uses `Simulation(...)` (flux-form convective step)
- LilyPad: uses `LilyPadSim(...)` (semi-Lagrangian predictor/corrector, fixed timestep)

Everything else (`sim_step!`, `sim_time`, field access, measurement utilities)
stays the same.

## Example: Flow Over a Circle (Same as WaterLily README)

This is the same circle setup from WaterLily's README so we can compare directly.

```julia
using WaterLily
using LilyPad

function circle_lp(n, m; Re=500, U=1)
    radius, center = m/8, m/2 - 1
    sdf(x, t) = sqrt(sum(abs2, x .- center)) - radius

    LilyPadSim((n, m),   # domain size
               (U, 0),   # domain velocity (& velocity scale)
               2radius;  # length scale
               ν=2radius*U/Re, body=AutoBody(sdf), Δt=2.5)
end

lp = circle_lp(3 * 2^7, 2^8)
sim_step!(lp, 10) # for example
```
The only difference is using the `LilyPadSim` constructor instead of `WaterLily.Simulation`.

### Comparison to WaterLily Results

The semi-Lagrangian method enables extremely large time steps without loss of stability, but this naturally changes details of the simulated solution. You can reproduce the README circle comparison locally by running:

```sh
julia --project=. examples/readme_circle_compare.jl gpu
```

Vorticity (`tU/L = 75`):

![Circle vorticity comparison](assets/circle_vorticity_compare.png)

Force histories:

![Circle force history comparison](assets/circle_force_compare.png)

Outputs of `examples/readme_circle_compare.jl` (`L=64`, `Re=500`), averaged over whole lift cycles (upward zero crossings) in `t=30–60`, timed on an RTX A2000 Laptop GPU.

|                | Cd    | Cd′   | CL    | St    | s/(tU/L) |
|----------------|-------|-------|-------|-------|----------|
| WaterLily      | 1.724 | 0.130 | 1.431 | 0.258 | 1.373    |
| LilyPad Δt=1.5 | 1.663 | 0.131 | 1.469 | 0.259 | 0.344    |
| LilyPad Δt=2.5 | 1.644 | 0.131 | 1.525 | 0.259 | 0.202    |

LilyPad is within 5% of WaterLily on Cd, 1% on Cd′ and St, and 7% on CL, and 4× faster at `Δt=1.5`, 7× faster at `Δt=2.5`. Larger steps lower the drag and raise the lift amplitude.

Each LilyPad step transports `uⁿ` along itself over `Δt` and projects (predictor), then transports `uⁿ - (Δt/2)∇p` along the average of `uⁿ` and the predicted velocity and projects again (corrector), so it costs two projections like WaterLily's step. Viscous diffusion acts along the trajectories, which needs `νΔt ≤ 1/3`; at this Reynolds number that caps the step near `Δt=2.6`.

## Quick Switching Pattern

If you already have a WaterLily setup function, make the constructor injectable:

```julia
using WaterLily
using LilyPad

function circle(make_sim, n, m; ν=0, U=1)
    radius, center = m/8, m/2 - 1
    sdf(x, t) = sqrt(sum(abs2, x .- center)) - radius
    make_sim((n, m), (U, 0), 2radius; ν, body=AutoBody(sdf))
end

wl = circle(Simulation, 96, 64)
lp = circle(LilyPadSim, 96, 64)
```

## Current Scope

This package is still under development.

- Validated in 2D and 3D on SIMD, CPU, and GPU.
- Viscous diffusion smooths each departure sample with a Gaussian of variance `2νΔt`; stable for `νΔt ≤ 1/3`.
- No acceleration term: `g` and `udf` are not supported.
- The time-stepping defaults at `Δt = 1.5`, but an adaptive scheme that takes advantage of SL-integration stability would be better.

For method details and development notes, see `CLAUDE.md`.

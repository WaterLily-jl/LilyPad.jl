module LilyPadBiotExt

using LilyPad
using WaterLily
using BiotSavartBCs

"""
    LilyBiotSim(dims, uBC, L; Δt=1.5, nonbiotfaces=(), fmm=true, symmetry=(), mem=Array, kwargs...)

Semi-Lagrangian simulation with Biot-Savart boundary conditions for unbounded external flow.
Combines `LilyFlow` (semi-Lagrangian advection) with `BiotSavartPoisson` (pressure solver).
Accepts the same arguments as `WaterLily.Simulation`.

The timestep `Δt` is the primary method knob: the semi-Lagrangian scheme has no
CFL constraint from advection, so larger values trade accuracy for speed.

- `nonbiotfaces`: tuple of face indices to exclude from Biot-Savart BCs (e.g. `(-2,)`)
- `fmm`: use the Fast Multi-level Method (`true`, default) or tree-sum (`false`)
- `symmetry`: tuple of symmetry-plane face indices (e.g. `(-2,)`)
"""
function LilyPad.LilyBiotSim(dims, uBC, L; Δt=1.5, nonbiotfaces=(), fmm=true, symmetry=(), mem=Array, kwargs...)
    WaterLily.Simulation(dims, uBC, L; flow_ctor=LilyFlow,
        pois_ctor=flow -> BiotSavartPoisson(flow; nonbiotfaces, fmm, symmetry, mem),
        Δt, mem, kwargs...)
end

end # module LilyPadBiotExt

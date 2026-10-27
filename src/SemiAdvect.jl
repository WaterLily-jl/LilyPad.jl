import WaterLily: @loop, loc
using Base: front

"""
    departure(x, u, dt, perdir=())

Characteristic backtrack for the departure point using RK2
"""
@inline function departure(x, u, dt, perdir=())
    xm = x - dt/2 * vinterp(x, u, perdir)
    return x - dt * vinterp(xm, u, perdir)
end

"""
    advect!(f, src, u, dt, V, perdir=(), ν=0)

Semi-Lagrangian advection of the vector field `src` by `u`, minus `V`, overwriting `f` in one kernel.
The departure sample is diffused over `dt` (Gaussian variance `2νdt`).
Faces outside `WaterLily.inside_u` (plus face 2 if periodic) are copied from `src`.
"""
advect!(f, src, u, dt, V, perdir=(), ν=0) = (σ² = 2ν * dt; @loop f[Ii] = _advect(Ii, src, u, dt, perdir, σ²) - V[Ii] over Ii ∈ CartesianIndices(f))
@inline function _advect(Ii::CartesianIndex{Np1}, src, u, dt, perdir, σ²) where Np1
    I, i, N = front(Ii), last(Ii), size(src)
    all(ntuple(d -> (d == i && d ∉ perdir ? 3 : 2) ≤ I[d] ≤ N[d] - (d == i), Np1 - 1)) || return src[Ii]
    quad(departure(loc(i, I, eltype(src)), u, dt, perdir), view(src, ntuple(_ -> :, Np1 - 1)..., i), i, perdir, σ²)
end

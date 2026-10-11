import WaterLily: @loop, loc

"""
    departure(x, u, dt, perdir=())

Characteristic backtrack for the departure point using RK2
"""
@inline function departure(x, u, dt, perdir=())
    xm = x - dt/2 * vinterp(x, u, perdir)
    return x - dt * vinterp(xm, u, perdir)
end

"""
    scalar_advect!(f, φ_src, u, dt, i, perdir=())

Semi-Lagrangian advection of face-component scalar field `φ_src` at face `i`, overwriting `f`.
"""
scalar_advect!(f::AbstractArray{T}, φ_src, u, dt, i, perdir=()) where T =
    @loop f[I] = quad(departure(loc(i, I, T), u, dt, perdir), φ_src, i, perdir) over I ∈ advect_range(size(f), i, perdir)
# `WaterLily.inside_u`, plus the face at index 2 if periodic
advect_range(N::NTuple{D}, i, perdir) where D = CartesianIndices(ntuple(d -> d == i ? ((i ∈ perdir ? 2 : 3):N[d]-1) : (2:N[d]), D))

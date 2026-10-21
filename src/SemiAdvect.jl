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
    advect!(f, (s, s⁻), (u, u⁻), dt, V, perdir=(), ν=0)

`f = (4s(x) - s⁻(x⁻))/3 - V`, where `x` backtracks one step under `u` and `x⁻` one more under `u⁻`,
and each sample is diffused over its age: `Δt` for `s(x)` and `2Δt` for `s⁻(x⁻)`.
"""
advect!(f, (s, s⁻), (u, u⁻), dt, V, perdir=(), ν=0) = @loop f[Ii] = _advect(Ii, s, s⁻, u, u⁻, dt, 2ν * dt, perdir) - V[Ii] over Ii ∈ CartesianIndices(f)
@inline function _advect(Ii::CartesianIndex{Np1}, s, s⁻, u, u⁻, dt, σ², perdir) where Np1
    I, i, N = front(Ii), last(Ii), size(s)
    all(ntuple(d -> (d == i && d ∉ perdir ? 3 : 2) ≤ I[d] ≤ N[d] - (d == i), Np1 - 1)) || return s[Ii]
    x = departure(loc(i, I, eltype(s)), u, dt, perdir)
    x⁻ = departure(x, u⁻, dt, perdir)
    (4quad(x, component(s, i), i, perdir, σ²) - quad(x⁻, component(s⁻, i), i, perdir, 2σ²)) / 3
end
@inline component(s::AbstractArray{T,Np1}, i) where {T,Np1} = view(s, ntuple(_ -> :, Np1 - 1)..., i)

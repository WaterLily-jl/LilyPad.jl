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
    advect!(f, s′, (s, s̃), u, dt, V, perdir=())

`f = (4s(x) - s̃(x))/3 - V` and `s′ = s(x)`, where `x` backtracks one step under `u`.
"""
advect!(f, s′, (s, s̃), u, dt, V, perdir=()) = @loop (f[Ii], s′[Ii]) = _advect(Ii, s, s̃, u, dt, V, perdir) over Ii ∈ CartesianIndices(f)
@inline function _advect(Ii::CartesianIndex{Np1}, s, s̃, u, dt, V, perdir) where Np1
    I, i, N = front(Ii), last(Ii), size(s)
    all(ntuple(d -> (d == i && d ∉ perdir ? 3 : 2) ≤ I[d] ≤ N[d] - (d == i), Np1 - 1)) || return (s[Ii] - V[Ii], s[Ii])
    x = departure(loc(i, I, eltype(s)), u, dt, perdir)
    A, B = quad(x, component(s, i), i, perdir), quad(x, component(s̃, i), i, perdir)
    ((4A - B) / 3 - V[Ii], A)
end
@inline component(s::AbstractArray{T,Np1}, i) where {T,Np1} = view(s, ntuple(_ -> :, Np1 - 1)..., i)

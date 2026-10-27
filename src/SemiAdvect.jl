import WaterLily: @loop, loc, δ
using Base: front

"""
    departure(x, u, dt, perdir=())

Characteristic backtrack for the departure point using RK2
"""
@inline function departure(x, u, dt, perdir=())
    xm = x - dt/2 * vinterp(x, u, perdir)
    return x - dt * vquad(xm, u, perdir)
end

"""
    advect!(f, s, u, dt, V, perdir=(), ν=0; disp)

`f = s(x) - V`, where `x` backtracks over `dt` under `u` and the sample is diffused over `dt`.
"""
function advect!(f::AbstractArray{T,Np1}, s, u, dt, V, perdir=(), ν=0;
                 disp=similar(f, SVector{Np1-1,T}, front(size(f)))) where {T,Np1}
    @loop disp[I] = _displacement(I, u, dt, perdir) over I ∈ CartesianIndices(disp)
    @loop f[Ii] = _advect(Ii, s, disp, 2ν * dt, perdir) - V[Ii] over Ii ∈ CartesianIndices(f)
end
@inline _displacement(I, u, dt, perdir) = (x = loc(0, I, eltype(u)); x - departure(x, u, dt, perdir))
@inline function _advect(Ii::CartesianIndex{Np1}, s, disp, σ², perdir) where Np1
    I, i, N = front(Ii), last(Ii), size(s)
    all(ntuple(d -> (d == i && d ∉ perdir ? 3 : 2) ≤ I[d] ≤ N[d] - (d == i), Np1 - 1)) || return s[Ii]
    quad(loc(i, I, eltype(s)) - toface(disp, I, i), s, i, perdir, σ², i)
end

@inline function toface(c, I, i)
    j = δ(i, I)
    2 ≤ I[i] - 2 && I[i] + 1 ≤ size(c, i) - 1 ? (9(c[I-j] + c[I]) - c[I-2j] - c[I+j]) / 16 : (c[I-j] + c[I]) / 2
end

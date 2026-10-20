import WaterLily
using WaterLily: _interp

# Clamp `x` to [lo,hi], or wrap it into [0,n) if periodic
@inline function _bound(x::T, lo, hi, n, per) where T
    per || return clamp(x, T(lo), T(hi))
    y = mod(x, T(n)); ifelse(y < n, y, zero(T))
end

# `WaterLily.interp` of vector field `u`, wrapping in the `perdir` directions
@inline vinterp(x::SVector{D,T}, u, perdir=()) where {D,T} = SVector{D,T}(_interp(SVector{D,T}(
    _bound(x[d] + T(d == i) / 2, 0, size(u, d) - 2, size(u, d) - 2, d ∈ perdir) for d in 1:D), view(u, ntuple(_ -> :, D)..., i)) for i in 1:D)

"""
    quad(x, f, i=0, perdir=())

Quadratic interpolation of the face-`i` field `f` at `x` (`i=0` for cell-centred `f`).
Queries are clamped to the domain, or wrapped in the `perdir` directions.
"""
@inline function quad(x::SVector{D,T}, f, i=0, perdir=()) where {D,T}
    xh = SVector{D,T}(_bound(x[d] + T(d == i) / 2, 0.5, size(f, d) - 2.5 + T(d == i), size(f, d) - 2, d ∈ perdir) for d in 1:D) .+ T(1.5)
    c = ntuple(d -> clamp(round(Int, xh[d]), 2, size(f, d) - 1), D)
    w = ntuple(d -> (t = xh[d] - c[d]; (t * (t - 1) / 2, 1 - t * t, t * (t + 1) / 2)), D)
    acc = zero(eltype(f))
    for J in CartesianIndices(ntuple(_ -> 3, D))
        @inbounds acc += f[CartesianIndex(ntuple(d -> c[d] + J.I[d] - 2, D))] * prod(ntuple(d -> w[d][J.I[d]], D))
    end
    return acc
end

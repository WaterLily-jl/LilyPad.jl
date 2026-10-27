import WaterLily
using WaterLily: _interp

# Wrap `x` into [0,n)
@inline _wrap(x::T, n) where T = (y = mod(x, T(n)); ifelse(y < n, y, zero(T)))
# Clamp `x` to [lo,hi], or wrap it into [0,n) if periodic
@inline _bound(x::T, lo, hi, n, per) where T = per ? _wrap(x, n) : clamp(x, T(lo), T(hi))

# `WaterLily.interp` of vector field `u`, wrapping in the `perdir` directions
@inline vinterp(x::SVector{D,T}, u, perdir=()) where {D,T} = SVector{D,T}(_interp(SVector{D,T}(
    _bound(x[d] + T(d == i) / 2, 0, size(u, d) - 2, size(u, d) - 2, d ∈ perdir) for d in 1:D), u, i) for i in 1:D)

@inline vquad(x::SVector{D,T}, u, perdir=()) where {D,T} = SVector{D,T}(quad(x, u, i, perdir, zero(T), i) for i in 1:D)

"""
    quad(x, f, i=0, perdir=(), σ²=0, c...)

Quadratic interpolation of the face-`i` field `f[..., c...]` at `x` (`i=0` for cell-centred `f`),
smoothed by a Gaussian of variance `σ²`. Beyond the ghosts, queries extrapolate linearly in the tangential directions and hold the last value in the normal direction `i`, or wrap in the `perdir` directions.
"""
@inline function quad(x::SVector{D,T}, f, i=0, perdir=(), σ²=zero(T), c...) where {D,T}
    xh = SVector{D,T}((y = x[d] + T(d == i) / 2; d ∈ perdir ? _wrap(y, size(f, d) - 2) : y) for d in 1:D) .+ T(1.5)
    I = ntuple(d -> clamp(round(Int, xh[d]), 2, size(f, d) - 1), D)
    w = ntuple(d -> (t = xh[d] - I[d]; t = d == i ? clamp(t, -one(T), one(T)) : t; (w₋, w₀, w₊) = t < -1 ? (-t, 1 + t, zero(T)) : t > 1 ? (zero(T), 1 - t, t) :
                     (t * (t - 1) / 2, 1 - t * t, t * (t + 1) / 2); (w₋ + σ² / 2, w₀ - σ², w₊ + σ² / 2)), D)
    acc = zero(eltype(f))
    for J in CartesianIndices(ntuple(_ -> 3, D))
        @inbounds acc += f[CartesianIndex(ntuple(d -> I[d] + J[d] - 2, D)), c...] * prod(ntuple(d -> w[d][J[d]], D))
    end
    return acc
end

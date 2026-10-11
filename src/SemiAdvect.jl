import WaterLily: @loop, inside_u, loc, interp

"""
    departure(x, u, dt)

Characteristic backtrack for the departure point using RK2
"""
@inline function departure(x, u, dt)
    xm = x - dt/2 * interp(x, u)
    return x - dt * interp(xm, u)
end

"""
    scalar_advect!(f, φ_src, u, dt, i)

Semi-Lagrangian advection of face-component scalar field `φ_src` at face `i`, overwriting `f`.
"""
scalar_advect!(f::AbstractArray{T}, φ_src, u, dt, i) where T = @loop f[I] = quad(departure(loc(i, I, T), u, dt), φ_src, i) over I ∈ inside_u(size(f), i)

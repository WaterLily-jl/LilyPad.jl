import WaterLily: @loop, inside_u, ∂, μddn, perBC!

"""
    LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}

Thin wrapper around `WaterLily.Flow` that overrides `mom_predict!`, `mom_correct!`
`BDIM!` and `CFL` with semi-Lagrangian advection functions.
Inviscid with no acceleration term: `ν≠0` warns, `g` and `udf` throw.
"""
struct LilyFlow{D, T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Tf<:AbstractArray{T}, Lf} <: AbstractFlow{D,T}
    flow :: Flow{D,T,Sf,Vf,Tf,Lf}
end
function LilyFlow(args...; kwargs...)
    flow = Flow(args...; kwargs...)
    isnothing(flow.g) || throw(ArgumentError("LilyFlow does not support `g`"))
    iszero(flow.ν) || @warn "LilyFlow is inviscid: ν=$(flow.ν) is ignored"
    LilyFlow(flow)
end

Base.getproperty(a::LilyFlow, name::Symbol) =
    name === :flow ? getfield(a, :flow) : getproperty(getfield(a, :flow), name)
Base.setproperty!(a::LilyFlow, name::Symbol, x) =
    setproperty!(getfield(a, :flow), name, x)

"""    mom_predict!(a::LilyFlow, t₀, t₁)

Semi-Lagrangian predictor: self-advect `u⁰` and apply BDIM blending.
"""
function WaterLily.mom_predict!(a::LilyFlow{D}, t₀, t₁; udf=nothing, kwargs...) where D
    isnothing(udf) || throw(ArgumentError("LilyFlow does not support `udf`"))
    advect!(a.f, a.u⁰, a.u⁰, t₁ - t₀, a.V, a.perdir)
    BDIM!(a)
end

"""    mom_correct!(a::LilyFlow, t)

Semi-Lagrangian corrector: advect `u⁰ - 0.5*dt*μ₀*∇p` under the averaged
velocity `0.5*(u⁰+u_pred)` and apply BDIM blending. Overwrites `u⁰` with the advected source.
"""
function WaterLily.mom_correct!(a::LilyFlow{D}, t; kwargs...) where D
    dt = a.Δt[end]
    a.u .= (a.u .+ a.u⁰) ./ 2
    @loop a.u⁰[Ii] -= (dt / 2) * a.μ₀[Ii] * ∂(last(Ii), front(Ii), a.p) over Ii ∈ inside_u(size(a.p))
    perBC!(a.u⁰, a.perdir)
    advect!(a.f, a.u⁰, a.u, dt, a.V, a.perdir)
    BDIM!(a)
end

"""
    CFL(a::LilyFlow)

Semi-Lagrangian step size is unlimited, so return the current timestep `a.Δt[end]`.
"""
WaterLily.CFL(a::LilyFlow; kwargs...) = a.Δt[end]

"""
    BDIM!(a::LilyFlow)

Semi-Lagrangian BDIM blending on the advected velocity field (`advect!` has already subtracted `V`).
"""
BDIM!(a::LilyFlow) = @loop a.u[Ii] = μddn(Ii, a.μ₁, a.f) + a.V[Ii] + a.μ₀[Ii] * a.f[Ii] over Ii ∈ inside_u(size(a.p))


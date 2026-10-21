import WaterLily: @loop, inside_u, μddn, mom_project!, CFL

"""
    LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}

`WaterLily.Flow` advanced by a semi-Lagrangian BDF2 step. Inviscid: `ν≠0` warns, `g` and `udf` throw.
The history `ũ` is `uⁿ⁻¹` carried forward to `tⁿ`, and `ũ′` is the buffer its update is written into.
"""
mutable struct LilyFlow{D, T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Tf<:AbstractArray{T}, Lf} <: AbstractFlow{D,T}
    const flow :: Flow{D,T,Sf,Vf,Tf,Lf}
    ũ :: Vf
    ũ′ :: Vf
end
function LilyFlow(args...; kwargs...)
    flow = Flow(args...; kwargs...)
    isnothing(flow.g) || throw(ArgumentError("LilyFlow does not support `g`"))
    iszero(flow.ν) || @warn "LilyFlow is inviscid: ν=$(flow.ν) is ignored"
    LilyFlow(flow, copy(flow.u⁰), similar(flow.u))
end

Base.getproperty(a::LilyFlow, name::Symbol) =
    name ∈ (:flow, :ũ, :ũ′) ? getfield(a, name) : getproperty(getfield(a, :flow), name)
Base.setproperty!(a::LilyFlow, name::Symbol, x) =
    name ∈ (:ũ, :ũ′) ? setfield!(a, name, x) : setproperty!(getfield(a, :flow), name, x)

"""
    mom_step!(a::LilyFlow, b::AbstractPoisson)

Semi-Lagrangian BDF2 step `(3uⁿ⁺¹ - 4uⁿ(x) + ũⁿ(x))/2Δt = -∇pⁿ⁺¹` with BDIM, where `x` backtracks
one step under `u⁺ = (3uⁿ - uⁿ⁻¹)/2` and `ũⁿ⁺¹ = uⁿ(x)` is kept as the next step's history.
"""
function WaterLily.mom_step!(a::LilyFlow{D,T}, b::WaterLily.AbstractPoisson; udf=nothing, kwargs...) where {D,T}
    isnothing(udf) || throw(ArgumentError("LilyFlow does not support `udf`"))
    @loop (a.u⁰[Ii], a.u[Ii]) = (a.u[Ii], (3a.u[Ii] - a.u⁰[Ii]) / 2) over Ii ∈ CartesianIndices(a.u)
    advect!(a.f, a.ũ′, (a.u⁰, a.ũ), a.u, a.Δt[end], a.V, a.perdir)
    BDIM!(a)
    a.p .*= T(2/3); mom_project!(a, b, 1, sum(a.Δt)); a.p .*= T(3/2)
    a.ũ, a.ũ′ = a.ũ′, a.ũ
    push!(a.Δt, CFL(a))
end

"""
    CFL(a::LilyFlow)

Fixed time step `a.Δt[end]`.
"""
WaterLily.CFL(a::LilyFlow; kwargs...) = a.Δt[end]

"""
    BDIM!(a::LilyFlow)

BDIM blending of the advected velocity `f`, which already has `V` subtracted.
"""
BDIM!(a::LilyFlow) = @loop a.u[Ii] = μddn(Ii, a.μ₁, a.f) + a.V[Ii] + a.μ₀[Ii] * a.f[Ii] over Ii ∈ inside_u(size(a.p))

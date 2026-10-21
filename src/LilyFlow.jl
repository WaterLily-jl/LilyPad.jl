import WaterLily: @loop, inside_u, μddn, mom_project!, CFL

"""
    LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}

`WaterLily.Flow` advanced by a semi-Lagrangian BDF2 step, stable for `νΔt ≤ 1/3`. `g` and `udf` throw.
"""
struct LilyFlow{D, T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Tf<:AbstractArray{T}, Lf} <: AbstractFlow{D,T}
    flow :: Flow{D,T,Sf,Vf,Tf,Lf}
    u⁺ :: Vf
    ū :: Vf
end
function LilyFlow(args...; kwargs...)
    flow = Flow(args...; kwargs...)
    isnothing(flow.g) || throw(ArgumentError("LilyFlow does not support `g`"))
    flow.ν * flow.Δt[end] ≤ 1/3 || throw(ArgumentError("LilyFlow needs νΔt ≤ 1/3, got $(flow.ν * flow.Δt[end])"))
    LilyFlow(flow, similar(flow.u), similar(flow.u))
end

Base.getproperty(a::LilyFlow, name::Symbol) =
    name ∈ (:flow, :u⁺, :ū) ? getfield(a, name) : getproperty(getfield(a, :flow), name)
Base.setproperty!(a::LilyFlow, name::Symbol, x) =
    setproperty!(getfield(a, :flow), name, x)

"""
    mom_step!(a::LilyFlow, b::AbstractPoisson)

Semi-Lagrangian BDF2 step `(3uⁿ⁺¹ - 4uⁿ(x) + uⁿ⁻¹(x⁻))/2Δt = -∇pⁿ⁺¹ + ν∇²u` with BDIM.
"""
function WaterLily.mom_step!(a::LilyFlow, b::WaterLily.AbstractPoisson; udf=nothing, kwargs...)
    isnothing(udf) || throw(ArgumentError("LilyFlow does not support `udf`"))
    @loop (a.u⁺[Ii], a.ū[Ii]) = ((3a.u[Ii] - a.u⁰[Ii]) / 2, (a.u[Ii] + a.u⁰[Ii]) / 2) over Ii ∈ CartesianIndices(a.u)
    advect!(a.f, (a.u, a.u⁰), (a.u⁺, a.ū), a.Δt[end], a.V, a.perdir, a.ν)
    a.u⁰ .= a.u
    BDIM!(a)
    a.p .*= eltype(a.p)(2/3) # projecting over Δt instead of 2Δt/3
    mom_project!(a, b, 1, sum(a.Δt))
    a.p .*= eltype(a.p)(3/2)
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

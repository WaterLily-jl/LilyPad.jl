import WaterLily: @loop, inside_u, μddn, mom_project!, CFL

"""
    LilyFlow{D,T,...} <: WaterLily.AbstractFlow{D,T}

Semi-Lagrangian advection-reflection `Flow`, stable for `νΔt ≤ 2/3`.
"""
struct LilyFlow{D, T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Tf<:AbstractArray{T}, Lf, Df<:AbstractArray{SVector{D,T}}} <: AbstractFlow{D,T}
    flow :: Flow{D,T,Sf,Vf,Tf,Lf}
    r :: Vf
    δp :: Sf
    disp :: Df
end
function LilyFlow(args...; kwargs...)
    flow = Flow(args...; kwargs...)
    isnothing(flow.g) || throw(ArgumentError("LilyFlow does not support `g`"))
    flow.ν * flow.Δt[end] ≤ 2/3 || throw(ArgumentError("LilyFlow needs νΔt ≤ 2/3, got $(flow.ν * flow.Δt[end])"))
    LilyFlow(flow, zero(flow.u), zero(flow.p), similar(flow.p, SVector{ndims(flow.p),eltype(flow.p)}))
end

Base.getproperty(a::LilyFlow, name::Symbol) =
    name ∈ (:flow, :r, :δp, :disp) ? getfield(a, name) : getproperty(getfield(a, :flow), name)
Base.setproperty!(a::LilyFlow, name::Symbol, x) = setproperty!(getfield(a, :flow), name, x)

"""
    mom_step!(a::LilyFlow, b::AbstractPoisson; udf=nothing)

Semi-Lagrangian advection-reflection step (Zehnder et al. 2018). Not compatible with a `udf` user-defined function.
"""
function WaterLily.mom_step!(a::LilyFlow, b::WaterLily.AbstractPoisson; udf=nothing, kwargs...)
    isnothing(udf) || throw(ArgumentError("LilyFlow does not support `udf`"))
    h, t = a.Δt[end] / 2, sum(a.Δt)
    @loop (a.u⁰[Ii], a.u[Ii]) = (a.u[Ii], (5a.u[Ii] - a.u⁰[Ii]) / 4) over Ii ∈ CartesianIndices(a.u)
    advect!(a.f, a.u⁰, a.u, h, a.V, a.perdir, a.ν; a.disp); BDIM!(a)
    a.p .-= a.δp; mom_project!(a, b, 2, t - h; itmx=1)
    r, u, u⁰, μ₀, p, perdir = a.r, a.u, a.u⁰, a.μ₀, a.p, a.perdir
    @loop (r[Ii], u[Ii]) = (reflect(Ii, u, μ₀, p, h, perdir), (3u[Ii] - u⁰[Ii]) / 2) over Ii ∈ CartesianIndices(u)
    advect!(a.f, a.r, a.u, h, a.V, a.perdir, a.ν; a.disp); BDIM!(a)
    correct!(a.u, a, 2h); a.p .+= a.δp; a.δp .= a.p .- a.δp
    mom_project!(a, b, 1, t); a.δp .= a.p .- a.δp
    push!(a.Δt, CFL(a))
end

@inline function reflect(Ii::CartesianIndex{Np1}, u, μ₀, p, h, perdir) where Np1
    I, i, N = front(Ii), last(Ii), size(p)
    i ∉ perdir && !(2 ≤ I[i] ≤ N[i] - 1) && return u[Ii]
    J = CartesianIndex(ntuple(d -> d ∈ perdir ? (I[d] == 1 ? N[d] - 1 : I[d] == N[d] ? 2 : I[d]) : clamp(I[d], 2, N[d] - 1), Np1 - 1))
    u[Ii] - h * μ₀[J, i] * WaterLily.∂(i, J, p)
end

function correct!(u, a::LilyFlow, c)
    μ₀, p = a.μ₀, a.p
    for i ∈ 1:ndims(p); @loop u[I, i] += c * μ₀[I, i] * WaterLily.∂(i, I, p) over I ∈ inside(p); end
end

"""
    CFL(a::LilyFlow)

Fixed time step `a.Δt[end]`.
"""
WaterLily.CFL(a::LilyFlow; kwargs...) = a.Δt[end]

"""
    BDIM!(a::LilyFlow)

BDIM blending of the advected velocity `f`, which has `V` subtracted.
"""
BDIM!(a::LilyFlow) = @loop a.u[Ii] = μddn(Ii, a.μ₁, a.f) + a.V[Ii] + a.μ₀[Ii] * a.f[Ii] over Ii ∈ inside_u(size(a.p))

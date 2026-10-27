using CUDA, WaterLily, LilyPad, Plots, DelimitedFiles, JLD2, Printf

function tgv(make_sim, p; Re=1600, T=Float32, kw...)
    L = 2^p; κ = T(π / L)
    function u0(i, xyz)
        x, y, z = @. xyz * κ
        i == 1 && return -sin(x) * cos(y) * cos(z)
        i == 2 && return cos(x) * sin(y) * cos(z)
        return zero(x)
    end
    make_sim((L, L, L), (0, 0, 0), 1 / κ; U=1, u0, ν=1 / (κ * Re), T, mem=CuArray, kw...)
end

function energy_enstrophy(sim)
    σ, u = sim.flow.σ, sim.flow.u; Ω = length(inside(σ))
    @inside σ[I] = WaterLily.ke(I, u)
    E = sum(@view σ[inside(σ)]) / Ω
    @inside σ[I] = WaterLily.ω_mag(I, u)
    E, sum(abs2, @view σ[inside(σ)]) * sim.flow.ν * sim.L / Ω
end

function history!(sim; t_max=20)
    rows = [(0.0, energy_enstrophy(sim)...)]
    while sim_time(sim) < t_max
        sim_step!(sim; remeasure=false)
        push!(rows, (sim_time(sim), energy_enstrophy(sim)...))
    end
    t, E, Z = (getindex.(rows, k) for k in 1:3)
    (; t, E, Z)
end

function cost(sim; t_max=20)
    sim_step!(sim; remeasure=false); CUDA.synchronize(); t₀ = sim_time(sim)
    wall = @elapsed (sim_step!(sim, t_max; remeasure=false); CUDA.synchronize())
    wall / (sim_time(sim) - t₀)
end

datadir = get(ENV, "TGV_DATA", joinpath(@__DIR__, "data", "tgv"))
dns = readdlm(joinpath(datadir, "TGV_Re1600.dat"), skipstart=43)
wl = jldopen(joinpath(datadir, "p7.jld2"))
pE = plot(dns[:, 1], dns[:, 2]; label="DNS", c=:black, xlabel="tπ/L", ylabel="kinetic energy", ylims=(0, 0.15))
pZ = plot(dns[:, 1], dns[:, 4]; label="DNS", c=:black, xlabel="tπ/L", ylabel="dissipation ν⟨ω²⟩", ylims=(0, 0.015))
plot!(pE, wl["t"], wl["E"]; label="WaterLily", c=1, ls=:dash)
plot!(pZ, wl["t"], wl["Z"]; label="WaterLily", c=1, ls=:dash)
row(name, t, E, Z, s...) = println(@sprintf("%-15s max ν⟨ω²⟩=%.4f at t=%.1f  E(20)=%.4f", name, maximum(Z), t[argmax(Z)], E[end]), s...)

s_wl = cost(tgv(Simulation, 7))
row("DNS", dns[:, 1], dns[:, 2], dns[:, 4]); row("WaterLily", wl["t"], wl["E"], wl["Z"], @sprintf("  %.3f s/(tU/L)", s_wl))
for (k, Δt) in enumerate((1f0, 2.5f0, 5f0, 7.5f0))
    name = "LilyPad Δt=$Δt"; r = history!(tgv(LilyPadSim, 7; Δt)); s = cost(tgv(LilyPadSim, 7; Δt))
    row(name, r.t, r.E, r.Z, @sprintf("  %.3f s/(tU/L) (%.1f×)", s, s_wl / s))
    plot!(pE, r.t, r.E; label=name, c=k + 1)
    plot!(pZ, r.t, r.Z; label=name, c=k + 1)
end
savefig(plot(pE, pZ; layout=(1, 2), size=(1300, 500), plot_title="TGV Re=1600, 128³"), joinpath(@__DIR__, "..", "assets", "tgv_p7.png"))

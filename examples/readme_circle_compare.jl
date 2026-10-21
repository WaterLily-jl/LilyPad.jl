using WaterLily, LilyPad, Plots, Statistics, Printf
gpu = "gpu" ∈ ARGS
gpu && using CUDA
mem = gpu ? CUDA.CuArray : Array
sync() = gpu ? CUDA.synchronize() : nothing

function circle(make_sim; L=32, Re=500, U=1, kw...)
    m = 4L; r, c = L / 2, m / 2 - 1
    make_sim((3m ÷ 2, m), (U, 0), L; ν=U * L / Re, body=AutoBody((x, t) -> √sum(abs2, x .- c) - r), T=Float32, mem, kw...)
end

function forces!(sim; t_end=60, dt=0.1)
    ts = 1:dt:t_end
    f = map(ts) do t
        sim_step!(sim, t; remeasure=false)
        WaterLily.total_force(sim) ./ (0.5sim.L * sim.U^2)
    end
    ts, first.(f), last.(f)
end

function metrics(ts, drag, lift; t₀=30, t₁=60)
    peaks = filter(k -> lift[k-1] < lift[k] ≥ lift[k+1] && t₀ ≤ ts[k] ≤ t₁, 2:length(lift)-1)
    R = first(peaks):last(peaks)-1
    (Cd=abs(mean(drag[R])), Cd′=std(drag[R]), CL=√2std(lift[R]), St=(length(peaks) - 1) / (ts[last(peaks)] - ts[first(peaks)]))
end

function cost(sim; window=10)
    t₀ = sim_time(sim); sync()
    wall = @elapsed (sim_step!(sim, t₀ + window; remeasure=false); sync())
    wall / (sim_time(sim) - t₀)
end

function vorticity(sim, title)
    @inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
    flood(Array(sim.flow.σ)[inside(sim.flow.σ)]; clims=(-10, 10), title, border=:none, legend=false)
    body_plot!(sim)
end

sims = ("WaterLily" => circle(Simulation), "LilyPad" => circle(LilyPadSim; Δt=2.5))
device = gpu ? CUDA.name(CUDA.device()) : "CPU, $(Threads.nthreads()) threads"
println("Re=500, ν=", last(sims[1]).flow.ν, " (WL) ", last(sims[2]).flow.ν, " (LP), LP Δt=", last(sims[2]).flow.Δt[end], ", ", device)

hist = [name => forces!(sim) for (name, sim) in sims]
secs = [name => cost(sim) for (name, sim) in sims]

out = joinpath(@__DIR__, "..", "assets")
pf = plot(xlabel="tU/L", ylabel="force coefficient", size=(900, 400), legend=:bottomleft)
for (name, (ts, drag, lift)) in hist
    plot!(pf, ts, [drag lift]; label=["drag $name" "lift $name"], ls=name == "LilyPad" ? :dash : :solid)
end
savefig(pf, joinpath(out, "circle_force_compare.png"))
savefig(plot((vorticity(sim, "$name, t=$(round(sim_time(sim), digits=1))") for (name, sim) in sims)...; layout=(2, 1), size=(900, 760)),
        joinpath(out, "circle_vorticity_compare.png"))

println(@sprintf("%-10s %7s %7s %7s %7s %10s", "", "Cd", "Cd′", "CL", "St", "s/(tU/L)"))
for ((name, h), (_, s)) in zip(hist, secs)
    m = metrics(h...)
    println(@sprintf("%-10s %7.3f %7.3f %7.3f %7.3f %10.3f", name, m.Cd, m.Cd′, m.CL, m.St, s))
end
println(@sprintf("LilyPad speed-up: %.2f×", last(secs[1]) / last(secs[2])))

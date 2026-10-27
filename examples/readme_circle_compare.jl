using WaterLily, LilyPad, Plots, Statistics, Printf
gpu = "gpu" ∈ ARGS
gpu && using CUDA
mem = gpu ? CUDA.CuArray : Array
sync() = gpu ? CUDA.synchronize() : nothing

function circle(make_sim; L=64, Re=500, U=1, kw...)
    m = 4L; r, c = L / 2, m / 2 - 1
    make_sim((3m ÷ 2, m), (U, 0), L; ν=U * L / Re, body=AutoBody((x, t) -> √sum(abs2, x .- c) - r), T=Float32, mem, kw...)
end

function forces!(sim; t_end=60)
    ts, drag, lift = Float64[], Float64[], Float64[]
    while sim_time(sim) < t_end
        sim_step!(sim; remeasure=false); F = WaterLily.total_force(sim) ./ (0.5sim.L * sim.U^2)
        push!(ts, sim_time(sim)); push!(drag, F[1]); push!(lift, F[2])
    end
    ts, drag, lift
end

function metrics(ts, drag, lift; t₀=30, t₁=60)
    k = findall(t -> t₀ ≤ t ≤ t₁, ts); l = lift .- mean(lift[k])
    up = filter(i -> l[i] < 0 ≤ l[i+1], k[1:end-1])
    tc = [ts[i] - l[i] * (ts[i+1] - ts[i]) / (l[i+1] - l[i]) for i in up]
    R = first(up)+1:last(up)
    (Cd=abs(mean(drag[R])), Cd′=std(drag[R]), CL=√2std(lift[R]), St=(length(tc) - 1) / (tc[end] - tc[1]))
end

function cost(sim; window=5, reps=3)
    minimum(1:reps) do _
        t₀ = sim_time(sim); sync()
        wall = @elapsed (sim_step!(sim, t₀ + window; remeasure=false); sync())
        wall / (sim_time(sim) - t₀)
    end
end

function vorticity(sim, title)
    @inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
    flood(Array(sim.flow.σ)[inside(sim.flow.σ)]; clims=(-10, 10), title, border=:none, legend=false)
    body_plot!(sim)
end

sims = ("WaterLily" => circle(Simulation), "LilyPad Δt=2.5" => circle(LilyPadSim; Δt=2.5), "LilyPad Δt=5" => circle(LilyPadSim; Δt=5))
device = gpu ? CUDA.name(CUDA.device()) : "CPU, $(Threads.nthreads()) threads"
println("Re=500, ν=", join((last(s).flow.ν for s in sims), ", "), " (WL, LP Δt=2.5, LP Δt=5), ", device)

hist = [name => forces!(sim) for (name, sim) in sims]
secs = [name => cost(sim) for (name, sim) in sims]

out = joinpath(@__DIR__, "..", "assets")
pf = plot(xlabel="tU/L", ylabel="force coefficient", size=(900, 400), legend=:bottomleft)
for (name, (ts, drag, lift)) in hist
    k = ts .≥ 1
    plot!(pf, ts[k], [drag[k] lift[k]]; label=["drag $name" "lift $name"], ls=name == "WaterLily" ? :solid : :dash)
end
savefig(pf, joinpath(out, "circle_force_compare.png"))
savefig(plot((vorticity(sim, "$name, t=$(round(sim_time(sim), digits=1))") for (name, sim) in sims)...; layout=(3, 1), size=(900, 1140)),
        joinpath(out, "circle_vorticity_compare.png"))

println(@sprintf("%-15s %7s %7s %7s %7s %10s", "", "Cd", "Cd′", "CL", "St", "s/(tU/L)"))
for ((name, h), (_, s)) in zip(hist, secs)
    m = metrics(h...)
    println(@sprintf("%-15s %7.3f %7.3f %7.3f %7.3f %10.3f", name, m.Cd, m.Cd′, m.CL, m.St, s))
end
foreach(((name, s),) -> println(@sprintf("%s speed-up: %.1f×", name, last(secs[1]) / s)), secs[2:end])

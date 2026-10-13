using LilyPad
using WaterLily
using Test
using StaticArrays
using BiotSavartBCs
import Pkg

#=
Array backends (mirrors WaterLily): LILYPAD_BACKENDS=cpu|cuda|all (comma-separated, default "cpu").
GPU backends are opt-in and installed on demand (they are not test dependencies), e.g.
    LILYPAD_BACKENDS=all julia --project -e 'using Pkg; Pkg.test()'
=#
const LILYPAD_BACKENDS = filter(!isempty, strip.(split(lowercase(get(ENV, "LILYPAD_BACKENDS", "cpu")), ',')))
(isempty(LILYPAD_BACKENDS) || any(b -> b ∉ ("cpu", "cuda", "all"), LILYPAD_BACKENDS)) &&
    throw(ArgumentError("LILYPAD_BACKENDS must be a comma-separated list of cpu|cuda|all, got \"$(get(ENV, "LILYPAD_BACKENDS", ""))\""))
_cpu = any(b -> b in ("cpu", "all"), LILYPAD_BACKENDS)
_cuda = any(b -> b in ("cuda", "all"), LILYPAD_BACKENDS)

# A requested GPU backend that cannot be installed or is not functional is skipped with a warning
_cuda && try
    Pkg.add("CUDA"); @eval using CUDA
    CUDA.functional() || (@warn "Requested CUDA is not functional; skipping"; global _cuda = false)
catch e
    @warn "Requested CUDA could not be installed; skipping" exception=e; global _cuda = false
end

arrays = []
_cpu && push!(arrays, Array)
_cuda && push!(arrays, CUDA.CuArray)
isempty(arrays) && error("No functional backend available")
@show arrays

include("maintests.jl")

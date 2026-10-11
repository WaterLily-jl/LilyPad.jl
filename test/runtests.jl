using LilyPad
using WaterLily
using Test
using StaticArrays
using BiotSavartBCs

_cuda = try (@eval using CUDA; CUDA.functional()) catch _ false end

function setup_backends()
    arrays = [Array]
    _cuda && push!(arrays, CUDA.CuArray)
    @show arrays
    return arrays
end

arrays = setup_backends()
include("maintests.jl")

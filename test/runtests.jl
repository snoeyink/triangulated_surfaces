using Test

include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

include("bitset.jl")
inclulde("bitsetoriented.jl")
include("geometry.jl")
include("unionfind.jl")
include("conflicts.jl")
include("bdry_loop.jl")
include("backtrack.jl")

using Test

include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

include("bitset.jl")
include("geometry.jl")
include("conflicts.jl")
include("bdry_loop.jl")



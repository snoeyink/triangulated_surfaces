using Test

include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

include('bitsets128.jl')
include("bitset.jl")
include("bitsetoriented.jl")
include("geometry.jl")
include("unionfind.jl")
include("conflicts.jl")
include("bdry_loop.jl")
include("backtrack.jl")

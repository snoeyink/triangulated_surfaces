using Test

include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

include("bitsets128.jl")
include("geometry.jl")
include("conflicts.jl")
include("packed_union_find.jl")


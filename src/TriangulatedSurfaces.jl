module TriangulatedSurfaces
using Random
using StaticArrays
using Base: @propagate_inbounds

export Point3D, Triangle, Surface, 
    tetrahedron_with_origin, add_random_interior_points!,
    has_colinear_triple, has_coplanar_quad,
    BitSet128, BitSetOriented128,
    precompute_conflicts, edge_index, ueindex, u_edge, triangle_index,
    enumerate_triangulated_surfaces

const D = true # debug flag
const MIN_VERTICES = 5
const MAX_VERTICES = 16
const N = 6


include("BitSets128.jl")
include("Points.jl")
include("OrientedET.jl")
include("PackedUnionFind.jl")
include("Backtrack.jl")


points = tetrahedron_with_origin(scale=4)
push!(points, Point3D(2, 1, 0))
count, surf = enumerate_triangulated_surfaces(points)
print(count, " triangulated surfaces with vertices ", points, "\n")
end # module

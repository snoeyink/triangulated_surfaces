module TriangulatedSurfaces
using Random
using StaticArrays
using Base: @propagate_inbounds

export Point3D, Triangle, Surface, 
       tetrahedron_with_origin, add_random_interior_points!,
    has_colinear_triple, has_coplanar_quad,
    BitSet128

const debugprint = true
const MIN_VERTICES = 5
const MAX_VERTICES = 16
const N = 6


include("BitSet128.jl")
include("Points.jl")
include("EdgesTriangles.jl")
include("BdryLoop.jl")
include("Backtrack.jl")


# When we choose a max triangle to start with, we can immediately prune all triangles that conflict with it and all triangles with higher indices.
# Note that we can search different max triangles in parallel. 



# --- Triangulated Surface Enumeration via Backtracking ---

"""
    enumerate_triangulated_surfaces(points::Vector{Point3D})

Enumerate all valid triangulated spheres on `points` using backtracking.
Returns a vector of `Surface` instances satisfying:
1. Manifold: every edge in 0 or 2 triangles
2. No intersection: edges don't intersect triangles except at shared vertices
3. Connected: path exists between any two vertices
4. Link: triangles around each vertex form a single loop

"""
function enumerate_triangulated_surfaces(points::Vector{Point3D})
    n = length(points)
    !(MIN_VERTICES <= n <= MAX_VERTICES) && throw(ArgumentError("Need $(MIN_VERTICES)<= n <= $(MAX_VERTICES) points"))
    
    # Precompute conflict matrix and mappings
    triangle_map, edgesets = precompute_conflicts(points)
    
    # Initialize per-vertex triangle counters and compute `min_idx`:
    # scan triangles from highest to lowest index, counting how many times
    # each vertex appears; when all n vertices have count >= 3, record the
    # triangle index `min_idx` where the last vertex crossed the threshold.
    v_tri_cnt = zeros(Int, n)
    v3_count = 0
    min_tri = 3n-6
    max_tri = length(triangle_map)
    v_count_tri(v) = (v_tri_cnt[v] += 1) == 3 ? 1 : 0
    @inbounds for i in 1:max_tri
        a, b, c = triangle_map[i]
        v3_count += v_count_tri(a)
        v3_count += v_count_tri(b) 
        v3_count += v_count_tri(c)
        if v3_count == n # all vertices have count >= 3
            min_tri = i
            break
        end
    end

    for t::Int16 in max_tri:-1:min_tri
        tmax, tmap, esets, tri_table = build_tri_table(n, t, triangle_map, edgesets)
        b = BdryLoop(n)
        init_loop!(b, tmax, tmap, esets)
        s = BdryStack()
        backtrack!(b, s, tri_table, esets, out)
    return results
end
end


end # module


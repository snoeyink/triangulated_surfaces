module TriangulatedSurfaces
using Random
using Base: @propagate_inbounds

export Point3D, Triangle, Surface, 
       tetrahedron_with_origin, add_random_interior_points!,
    has_colinear_triple, has_coplanar_quad,
    BitSet128

const debugprint = true


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
    !(4<= n <= 16) && throw(ArgumentError("Need 4<= n <= 16 points"))
    
    # Precompute conflict matrix and mappings
    triangle_map, edgesets = precompute_conflicts(points)
     
    # Initialize per-vertex triangle counters and compute `min_idx`:
    # scan triangles from highest to lowest index, counting how many times
    # each vertex appears; when all n vertices have count >= 3, record the
    # triangle index `min_idx` where the last vertex crossed the threshold.
    v_tri_cnt = zeros(Int, n)
    v3_count = 0
    v_count_tri(v) = (v_tri_cnt[v] += 1) == 3 ? 1 : 0
    @inbounds for i in 1:length(tri_map)
        a, b, c = tri_map[i]
        v3_count += v_count_tri(a)
        v3_count += v_count_tri(b) 
        v3_count += v_count_tri(c)
        if v3_count == n
            min_tri = i
            break
        end
    end

    max_tri = triangle_index(n-2, n-1, n) # all triangles must be below this
    for t in max_tri:-1:min_tri
        tmax, tmap, esets, tri_table = build_tri_table(n, t, triangle_map, edgesets)
        results = backtrack(tmax, tmap, esets, tri_table)
    return results
end
end


end # module


module TriangulatedSurfaces
using Random
using Base: @propagate_inbounds

export Point3D, Triangle, Surface, 
       tetrahedron_with_origin, add_random_interior_points!,
    has_colinear_triple, has_coplanar_quad,
    BitSet128

include("BitSet128.jl")
include("Points.jl")
include("EdgesTriangles.jl")
include("BdryLoop.jl")
include("Backtrack.jl")


# When we choose a max triangle to start with, we can immediately prune all triangles that conflict with it and all triangles with higher indices.
# Note that we can search different max triangles in parallel. 
function initialize_with_maxtriangle(n::Int, max_tri::Int, triangle_map, ETConflict, edgeset)
    max_edge_idx = edge_index(n - 1, n)
    max_tri_idx = triangle_index(n - 2, n - 1, n)
    
    (a,b,c) = triangle_map[max_tri] # note: a<b<c always, so no need to check order for edge_bit
    curr_edgeset = singleton(a,b)|singleton(b,c)|singleton(a,c) 
    
    # reorder triidx by increasing conflict count for better pruning
    function triangle_conflict(mt, t) 
        edgesets_disjoint(curr_edgeset, edgeset[t]) || return true
        (a,b,c) = triangle_map[t]
        ETConflict[edge_index(a,b), mt] && return true
        ETConflict[edge_index(b,c), mt] && return true 
        ETConflict[edge_index(a,c), mt] && return true
        return false
    end
    
    tri_indices = [t for t in 1:max_tri_idx if ((t <= max_tri) && !triangle_conflict(max_tri, t))]
    t_index_map = zeros(Int, 16, 16, 16)
    @inbounds for (i, tri_idx) in enumerate(tri_indices)
        c, d, e = triangle_map[tri_idx]
        @inbounds for (x,y,z) in ((c,d,e),(c,e,d),(d,c,e),(d,e,c),(e,c,d),(e,d,c))
            t_index_map[x,y,z] = i
        end
    end
    # Build conflict lists from ETConflict matrix
    conflicted_edges(tri_idx) = [edge_idx for edge_idx in 1:max_edge_idx if ETConflict[edge_idx, tri_idx]]
    te_conflicts = [conflicted_edges(tri_idx) for tri_idx in tri_indices]

    return curr_edgeset, te_conflicts, t_index_map
end


function set_max_triangle(max_tri::Int, triangle_map, ETConflict)
    (a,b,c) = triangle_map[max_tri] # note: a<b<c always, so no need to check order for edge_bit
    curr_edgeset = edge_bit(a,b)|edge_bit(b,c)|edge_bit(a,c) 
    tri_indices = [t for t in 1:length(triangle_map) if ((t <= max_tri) && !triangle_conflict(max_tri, t))]
    return curr_edgeset, tri_indices
end

# --- Triangulated Surface Enumeration via Backtracking ---

"""
    enumerate_triangulated_surfaces(points::Vector{Point3D})

Enumerate all valid triangulated spheres on `points` using backtracking.
Returns a vector of `Surface` instances satisfying:
1. Manifold: every edge in 0 or 2 triangles
2. No intersection: edges don't intersect triangles except at shared vertices
3. Connected: path exists between any two vertices
4. Link: triangles around each vertex form a single loop

Limited to `max_results` surfaces to prevent exhaustive search.
"""
function enumerate_triangulated_surfaces(points::Vector{Point3D})
    n = length(points)
    !(4<= n <= 16) && throw(ArgumentError("Need 4<= n <= 16 points"))
    
    # Precompute conflict matrix and mappings
     conflictcount, edgeset, triangle_map, ETConflict = precompute_conflicts(points)
     
    # Initialize per-vertex triangle counters and compute `min_idx`:
    # scan triangles from highest to lowest index, counting how many times
    # each vertex appears; when all n vertices have count >= 3, record the
    # triangle index `min_idx` where the last vertex crossed the threshold.
    v_tri_cnt = zeros(Int, n)
    v3_count = 0
    v_count_tri(v) = (v_tri_cnt[v] += 1) == 3 ? 1 : 0
    @inbounds for i in length(tri_map):-1:1
        a, b, c = tri_map[i]
        v3_count += v_count_tri(a)
        v3_count += v_count_tri(b) 
        v3_count += v_count_tri(c)
        if v3_count == n
            min_tri = i
            break
        end
    end

    max_tri = triangle_index(1,2,n) # all triangles must be below this
    for t in min_tri:max_tri
        tri_map, ETConf, curr_edgeset, eset, te_conf, ti_map = initialize_with_maxtriangle(max_tri, conflictcount, edgeset, triangle_map, ETConflict, points)
    
    
    return results
end
end


end # module



"""
    conflict(edge, tri, points)

Return `true` when edge `(a,b)` intersects triangle `(c,d,e)` anywhere other than
sharing vertices. Uses an integer-only segment–triangle intersection test.
"""
function conflict(edge::Tuple{Int,Int}, tri::NTuple{3,Int}, points::Vector{Point3D})
    a, b = edge
    c, d, e = tri
    if a in tri || b in tri
        return false
    end
    pa, pb = points[a], points[b]
    pc, pd, pe = points[c], points[d], points[e]
    return segment_intersects_triangle(pa, pb, pc, pd, pe)
end

"""
    te_index(a::Int, b::Int, c::Int)
"""
triangle_index(a::Int, b::Int, c::Int) = a + (b-1)*(b-2) ÷ 2 + (c-1)*(c-2)*(c-3) ÷ 6 # 1 based

@inline function edge_index(a::Int, b::Int)
    @boundscheck (1 <= a < b <= 16) || throw(ArgumentError("edge endpoints must satisfy 1 <= a < b <= 16"))
    return a + (b-1)*(b-2) ÷ 2 # 1 based
end
e_index(a::Int, b::Int) = @inbounds edge_index(minmax(a,b)...) # 1 based
singleton(a::Int, b::Int) = @inbounds singleton(edge_index(a,b))

"""
    precompute_conflicts(points::Vector{Point3D})

Precompute conflicts:
1. triangle_map: map from triangle numbers 1:n choose 3 to triples of vertex indices (c,d,e) with c<d<e
1. ETConflict: Boolean matrix where ETConflict[edge_index(a,b), triangle_index(c,d,e)] 
    is true iff edge (a,b) conflicts with triangle (c,d,e)
"""
function precompute_conflicts(points::Vector{Point3D})
    n = length(points)
    max_edge_idx = edge_index(n - 1, n)
    max_tri_idx = triangle_index(n - 2, n - 1, n)
    
    # Create reverse map: index -> triangle
    triangle_map = Vector{NTuple{3,Int}}(undef, max_tri_idx)
    for c in 1:(n-2), d in (c+1):(n-1), e in (d+1):n
        triangle_map[triangle_index(c, d, e)] = (c, d, e)
    end
    
    # Initialize Boolean matrix for conflicts
    ETConflict = falses(max_edge_idx, max_tri_idx)
   
    # Populate conflict matrix: loop over edge and triangle numbers
    conflictcount = zeros(Int, max_tri_idx) 
    edgeset = fill(BitSet128(), max_tri_idx) # bitset of conflicting edges for each triangle
    for a in 1:(n-1), b in (a+1):n
        edge_idx = edge_index(a, b)
        edge_singleton = singleton(edge_idx)
        for tri_idx in 1:max_tri_idx
            conf = conflict((a, b), triangle_map[tri_idx], points) 
            ETConflict[edge_idx, tri_idx] = conf
            if conf
                conflictcount[tri_idx] += 1
                edgeset[tri_idx] |= edge_singleton
            end
        end
    end

    tri_indices = sortperm(1:max_tri_idx, by=i->conflictcount[i])
    print(conflictcount[tri_indices]) 
    return triangle_map[tri_indices], ETConflict[:, tri_indices], edgeset[tri_indices] 
end

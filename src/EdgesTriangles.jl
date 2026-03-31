
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


struct Tri_Edgesets
    has::BitSet128 # set of 3 edges in the triangle
    conf::BitSet128 # set of edges that conflict with the triangle
end

"""
    triangle_index(a::Int, b::Int, c::Int)
"""
triangle_index(a::Int, b::Int, c::Int) = a + (b-1)*(b-2) ÷ 2 + (c-1)*(c-2)*(c-3) ÷ 6 # 1 based

@inline function edge_index(a::Int, b::Int)
    @boundscheck (1 <= a < b <= 16) || throw(ArgumentError("edge endpoints must satisfy 1 <= a < b <= 16"))
    return a + (b-1)*(b-2) ÷ 2 # 1 based
end
e_index(a::Int, b::Int) = @inbounds edge_index(minmax(a,b)...) 
singleton(a::Int, b::Int) = @inbounds singleton(e_index(a,b))

"""
    precompute_conflicts(points::Vector{Point3D})

Precompute conflicts:
1. triangle_map: map from triangle numbers 1:n choose 3 to triples of vertex indices (c,d,e) with c<d<e
2. edgesets: parallel array with BitSet128s for edges triangle has and conflicts with. 
These are sorted by increasing conflicts 
"""
function precompute_conflicts(points::Vector{Point3D})
    n = length(points)
    max_edge_idx = edge_index(n - 1, n)
    max_tri_idx = triangle_index(n - 2, n - 1, n)
    
    # Create reverse map: index -> triangle
    triangle_map = Vector{NTuple{3,Int}}(undef, max_tri_idx)
    for c in 1:(n-2), d in (c+1):(n-1), e in (d+1):n
        t = triangle_index(c, d, e)
        triangle_map[t] = (c, d, e)
    end
    
    # Loop over triangles and edges to count conflicts and make edgesets
    conflictcount = zeros(Int, max_tri_idx) 
    edgesets = Vector{Tri_Edgesets}(undef, max_tri_idx) # bitset of conflicting edges for each triangle
    for t in 1:max_tri_idx
        conf_edgeset = BitSet128()
        for a in 1:(n-1), b in (a+1):n 
            if conflict((a, b), triangle_map[t], points)
                conflictcount[t] += 1
                conf_edgeset |= singleton(a,b)
                for c in 1:a-1 conflictcount[triangle_index(c, a, b)] += 1; end
                for c in a+1:b-1 conflictcount[triangle_index(a, c, b)] += 1; end
                for c in b+1:n conflictcount[triangle_index(a, b, c)] += 1; end
            end
        end # for edge
        v1, v2, v3 = triangle_map[t]
        edgesets[t] = Tri_Edgesets(singleton(v1, v2) | singleton(v2, v3) | singleton(v1, v3),conf_edgeset) 
    end # for tri

    tri_indices = sortperm(1:max_tri_idx, by=i->conflictcount[i])
    debugprint && println(conflictcount[tri_indices])
    return triangle_map[tri_indices], edgesets[tri_indices] 
end

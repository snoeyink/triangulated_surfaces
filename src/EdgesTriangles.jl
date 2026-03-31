
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


# ── Tri-table builder ─────────────────────────────────────────────────────────
# Given a tmax triangle index, compacts triangle_map and edgesets 
# to include only triangles that are compatible with tmax's edgesets, 
# and builds tri_table for the remaining triangles.
# Since no triangle conflicting with tmax remains, 
# I can set esets[tmax].conflicts = esets[tmax].has so it conflicts with added_edgeset and serves as sentinel. 
# Fills tri_table[a,b,c] (all 6 permutations) with the renumbered index for
# triangles ≤ tmax that are compatible with esets[tmax].
# Compaction helps tmap & esets fit into L1 cache. 

function build_tri_table(n::Int,       # number of points 
                         tmax ::Int,   # max triangle is first used in surface
                         triangle_map     ::Vector{NTuple{3,Int}},
                         edgesets         ::Vector{Tri_Edgesets})
    @inbounds begin
        has_tmax = edgesets[tmax].has
        conf_tmax = edgesets[tmax].conf
        keep = [
            isdisjoint(edgesets[t].conf, has_tmax) &&
            isdisjoint(edgesets[t].has, conf_tmax)
            for t in 1:tmax
        ]
        tmap = triangle_map[keep]
        esets = edgesets[keep]
        tmax = length(esets)
        esets[tmax] = Tri_Edgesets(has_tmax, has_tmax) # make tmax conflict with itself so it won't be used again
    end

    tri_table = Array{Int16,3}(undef, n, n, n)
    fill!(tri_table, Int16(tmax)) # tmax is sentinel, too, because it conflictsts own edgeset
    @inbounds for t in 1:tmax
        a, b, c = tmap[t]
        idx = Int16(t)
        tri_table[a,b,c]=idx; tri_table[a,c,b]=idx
        tri_table[b,a,c]=idx; tri_table[b,c,a]=idx
        tri_table[c,a,b]=idx; tri_table[c,b,a]=idx
    end

    return tmax, tmap, esets, tri_table
end
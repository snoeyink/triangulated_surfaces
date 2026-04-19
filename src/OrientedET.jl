# OrientedET.jl

TriIJK = NTuple{3,UInt8} # triangle map
EdgeIJ = NTuple{2,UInt8} # edge map

"""
    conflict(edge, tri, points)

Return `true` when edge `(a,b)` intersects triangle `(c,d,e)` anywhere other than
sharing vertices. Uses an integer-only segment–triangle intersection test.
"""
function conflict(edge::Tuple{Integer,Integer}, tri::TriIJK, points::Vector{Point3D})
    a, b = edge
    c, d, e = tri
    if a in tri || b in tri
        return (a==c && b==d) || (a==d && b==e) || (a==e && b==c) # conflict with contained oriented edges
    end
    pa, pb = points[a], points[b]
    pc, pd, pe = points[c], points[d], points[e]
    return segment_intersects_triangle(pa, pb, pc, pd, pe)
end

"""
    edge_index(a::Integer, b::Integer)
    
Returns a strictly 1-based index where paired oriented edges (a,b) and (b,a) 
occupy adjacent indices: (1, 2), (3, 4), etc.
"""
@inline function edge_index(a::Integer, b::Integer)
    @boundscheck (1 <= a <= N) && (1 <= b <= N) && (a ≠ b) ||
        throw(ArgumentError("edge endpoints must satisfy 1 <= a, b <= $(N)"))
    return (a < b) ? 2*a + (b-1)*(b-2) - 1 : 2*b + (a-1)*(a-2)
end

"""
    rev_edge(idx::Integer)

Zero-cost mapping from an edge index to its reverse edge index using XOR.
"""
@inline rev_edge(idx::Integer) = ((idx - 1) ⊻ 1) + 1

singleton(a::Integer, b::Integer) = @inbounds singleton(edge_index(a,b))

"""
    triangle_index(a::Integer, b::Integer, c::Integer)
"""
@inline function triangle_index(a::Integer, b::Integer, c::Integer)
    @boundscheck (1 <= a <= N) && (1 <= b <= N) && (a ≠ b) && (1 <= c <= N) && (a ≠ c) && (b ≠ c)||
        throw(ArgumentError("edge endpoints must satisfy 1 <= a < b <= $(N)"))
    _,i = findmin((a,b,c))
    i == 1 && return @inbounds (a-1)*N*(N-1)+edge_index(b,c)
    i == 2 && return @inbounds (b-1)*N*(N-1)+edge_index(c,a)
    return @inbounds (c-1)*N*(N-1)+edge_index(a,b)
end

"""
    precompute_conflicts(points::Vector{Point3D})

Precompute conflicts:
1. triangle_map: map from triangle numbers 1:P(n,3) to triples of vertex indices (c,d,e) with c<d<e
2. edge_map: map from edge indices 1:P(n,2) to pairs of vertices
3. EconflT: BitSet128s for edges triangle conflicts with. 
These are sorted by increasing conflicts 
"""
function precompute_conflicts(points::Vector{Point3D})
    length(points) == N || throw(ArgumentError("points length must be equal to N=$(N)"))
    max_edge_idx = edge_index(N, N - 1) # This naturally maxes out at N*(N-1)
    max_tri_idx = triangle_index(N - 2, N, N - 1)

    # Create reverse map: index -> triangle
    triangle_map = Vector{TriIJK}(undef, max_tri_idx)
    edge_map = Vector{EdgeIJ}(undef, max_edge_idx)
    for a in 1:N, b in 1:N
        a == b && continue
        edge_map[edge_index(a,b)] = (a, b)
    end

    for c in 1:(N-2), d in (c+1):(N-1)
        for e in (d+1):N
            t = triangle_index(c, d, e)
            triangle_map[t] = (c, d, e)
            triangle_map[t+N*(N-1)÷ 2] = (c, e, d)
        end
    end
    
    # Loop over triangles and edges to count conflicts and make edgesets
    conflictcount = zeros(Int, max_edge_idx) 
    edges_conflT = [Set{Int}() for _ in 1:max_tri_idx]
    
    for a in 1:N, b in 1:N 
        a == b && continue
        ei = edge_index(a,b)
        confl_count = 0
        for t in 1:max_tri_idx
            if conflict((a, b), triangle_map[t], points)
                confl_count += 1
                push!(edges_conflT[t], ei)
            end
        end # for triangle
        conflictcount[ei] = confl_count 
    end # for edge

    EconflT = Vector{BitSet128}(undef, max_tri_idx) # conflicted edges for each triangle index
    
    # Stable sort preserves the (2k-1, 2k) pairing exactly!
    edge_indices = sortperm(1:max_edge_idx, by=i->conflictcount[i])
    ip = invperm(edge_indices)
    
    for t = 1:max_tri_idx
        conf_edgeset = BitSet128()
        for i in edges_conflT[t]
            conf_edgeset |= singleton(BitSet128, ip[i])
        end
        EconflT[t] = conf_edgeset
    end
    return triangle_map, edge_map[edge_indices], EconflT 
end
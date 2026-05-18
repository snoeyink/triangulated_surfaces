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
    
Returns a strictly 1-based oriented index with orientation as the high block:
- forward edges `(a,b)` with `a<b` map to `1:128`
- reverse edges `(a,b)` with `a>b` map to `129:256`

For a fixed undirected edge `{u,v}` with `u<v`, `(u,v)` and `(v,u)` share the
same low-order index and differ by `+128`.
"""
@inline function edge_index(a::Integer, b::Integer)
    @boundscheck (1 <= a <= N) && (1 <= b <= N) && (a ≠ b) ||
        throw(ArgumentError("edge endpoints must satisfy 1 <= a, b <= $(N)"))
    return (a < b) ? a + ((b-1)*(b-2)) ÷ 2 : b + ((a-1)*(a-2)) ÷ 2 + 128
end

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
2. edge_map: map from oriented edge indices 1:(128 + P(n,2)) to pairs of vertices
3. EconflT: BitSet128 for unoriented edges each triangle conflicts with

No edge reindexing is applied.
"""
function precompute_conflicts(points::Vector{Point3D})
    length(points) == N || throw(ArgumentError("points length must be equal to N=$(N)"))
    max_edge_idx = edge_index(N, N - 1) # max used index = 128 + binomial(N,2)
    max_tri_idx = 0
    for a in 1:N, b in 1:N, c in 1:N
        (a == b || b == c || a == c) && continue
        t = triangle_index(a, b, c)
        t > max_tri_idx && (max_tri_idx = t)
    end

    # Create reverse map: index -> triangle
    triangle_map = fill((UInt8(0), UInt8(0), UInt8(0)), max_tri_idx)
    edge_map = Vector{EdgeIJ}(undef, max_edge_idx)
    for a in 1:N, b in 1:N
        a == b && continue
        edge_map[edge_index(a,b)] = (a, b)
    end

    for a in 1:N, b in 1:N, c in 1:N
        (a == b || b == c || a == c) && continue
        triangle_map[triangle_index(a, b, c)] = (a, b, c)
    end
    
    # Loop over triangles and edges to count conflicts and make edgesets
    edges_conflT = [Set{Int}() for _ in 1:max_tri_idx]
    
    for a in 1:N, b in 1:N 
        a == b && continue
        ei = edge_index(a,b)
        for t in 1:max_tri_idx
            triangle_map[t][1] == 0 && continue
            if conflict((a, b), triangle_map[t], points)
                push!(edges_conflT[t], (ei <= 128 ? ei : ei - 128))
            end
        end # for triangle
    end # for edge

    EconflT = Vector{BitSet128}(undef, max_tri_idx) # conflicted unoriented edges for each triangle index
    
    for t = 1:max_tri_idx
        conf_edgeset = BitSet128()
        for i in edges_conflT[t]
            conf_edgeset |= singleton(BitSet128, i)
        end
        EconflT[t] = conf_edgeset
    end
    return triangle_map, edge_map, EconflT 
end

@inline undirected_edge_index(e::Int) = e <= 128 ? e : e - 128
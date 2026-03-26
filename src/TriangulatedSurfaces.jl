module TriangulatedSurfaces
using Random

export Point3D, Triangle, Surface, 
       tetrahedron_with_origin, add_random_interior_points!,
    has_colinear_triple, has_coplanar_quad,
    BitSet128

include("BitSet128.jl")

"""
    Point3D(x, y, z)

Minimal 3D point type.
"""
struct Point3D
    x::Int
    y::Int
    z::Int
end

"""
    Triangle(a, b, c)

Triangle by vertex indices (1-based into the surface's point list).
"""
struct Triangle
    a::Int
    b::Int
    c::Int
end

"""
    Surface(points, triangles)

Collection of vertices and their triangulation.
"""
struct Surface
    points::Vector{Point3D}
    triangles::Vector{Triangle}
end

# Lightweight vector helpers to keep code simple; all integer arithmetic.
@inline sub(a::Point3D, b::Point3D) = (a.x - b.x, a.y - b.y, a.z - b.z)
@inline dotp(u::NTuple{3,<:Integer}, v::NTuple{3,<:Integer}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
@inline @inbounds crossp(u::NTuple{3,<:Integer}, v::NTuple{3,<:Integer}) = (
    u[2] * v[3] - u[3] * v[2],
    u[3] * v[1] - u[1] * v[3],
    u[1] * v[2] - u[2] * v[1]
)
@inline norm2(u::NTuple{3,<:Integer}) = dotp(u, u)



"""
    tetrahedron_with_origin(; scale=4096)

Return five points: four vertices of a regular tetrahedron scaled by `scale`
and the origin. The origin is strictly inside this tetrahedron.
"""
function tetrahedron_with_origin(; scale::Int=4096)
    v1 = Point3D(scale, scale, scale)
    v2 = Point3D(scale, -scale, -scale)
    v3 = Point3D(-scale, scale, -scale)
    v4 = Point3D(-scale, -scale, scale)
    o = Point3D(0, 0, 0)
    return [v1, v2, v3, v4, o]
end

"""
    has_colinear_triple(points)

Return `true` if any triple of points is colinear (exact integer arithmetic)
"""
function has_colinear_triple(points::Vector{Point3D})
    n = length(points)
    n < 3 && return false
    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        u = sub(points[j], points[i])
        v = sub(points[k], points[i])
        if norm2(crossp(u, v)) <= 0
            return true
        end
    end
    return false
end

"""
    has_coplanar_quad(points)

Return `true` if any 3 points with the last point is coplanar (exact integer arithmetic)
"""
function has_coplanar_quad(points::Vector{Point3D})
    n = length(points)
    n < 4 && return false
    for i in 1:(n - 3), j in (i + 1):(n - 2), k in (j + 1):(n - 1)
        u = sub(points[j], points[i])
        v = sub(points[k], points[i])
        w = sub(points[n], points[i])
        if dotp(crossp(u, v), w) == 0
            return true
        end
    end
    return false
end

"""
    point_in_tetra(p, v1, v2, v3, v4)

Return `true` if `p` lies strictly inside the tetrahedron defined by
vertices `v1`..`v4`. Points on the boundary are rejected.
"""
function point_in_tetra(p::Point3D, v1::Point3D, v2::Point3D, v3::Point3D, v4::Point3D)
    orient(a, b, c, d) = dotp(crossp(sub(b, a), sub(c, a)), sub(d, a))
    s = sign(orient(v1, v2, v3, v4))
    return s == sign(orient(v2, v4, v3, p)) && 
           s == sign(orient(v1, v3, v4, p)) && 
           s == sign(orient(v1, v4, v2, p)) && 
           s == sign(orient(v1, v2, v3, p)) 
end


"""
    random_interior_points(tetra_vertices, n; rng=Random.default_rng())

Select `n` distinct random integer points strictly inside the tetrahedron
using rejection sampling.
"""
function add_random_interior_points!(tetra::Vector{Point3D}, n::Int; rng::AbstractRNG=Random.default_rng(), scale::Int=4096, max_attempts::Int=1000)
    n < 13 || throw(ArgumentError("n must be less than 13"))
    v1, v2, v3, v4, o = tetra
    plane(v1, v2, v3) = crossp(sub(v2, v1), sub(v3, v1))
    orient(plane, p, v) = dotp(plane, sub(p, v))
    # planes for tetrahedral faces
    p123 = plane(v1, v2, v3)
    p243 = plane(v2, v4, v3)
    p134 = plane(v1, v3, v4)
    p142 = plane(v1, v4, v2)
    s = sign(orient(p123, v4, v1)) # orientation of tetrahedron
    
    result = tetra
    attempts = 0
    while length(result) < n && attempts < max_attempts
        attempts += 1
        x = rand(rng, -scale:scale)
        y = rand(rng, -scale:scale)
        z = rand(rng, -scale:scale)
        p = Point3D(x, y, z)
        s123 = orient(p123, p, v1)
        s243 = orient(p243, p, v2)
        s134 = orient(p134, p, v1)
        s142 = orient(p142, p, v1)
        
        if sign(orient(p123, p, v1)) != s || sign(orient(p243, p, v2)) != s ||
           sign(orient(p134, p, v1)) != s || sign(orient(p142, p, v1)) != s
            continue
        end
        
        if !has_coplanar_quad(vcat(result, p))
            push!(result, p)
        end
    end
    return result    
end


# Segment–triangle intersection using only integer signed volumes.
# Assumes non-degenerate input (no coplanar quadruples).
@inline signed_volume6(a::Point3D, b::Point3D, c::Point3D, d::Point3D) = dotp(crossp(sub(b, a), sub(c, a)), sub(d, a))

function segment_intersects_triangle(p0::Point3D, p1::Point3D, v0::Point3D, v1::Point3D, v2::Point3D)
    # Endpoints on opposite sides of the triangle plane
    s0 = signed_volume6(v0, v1, v2, p0)
    s1 = signed_volume6(v0, v1, v2, p1)
    if (s0 == 0 || s1 == 0) print("Warning: segment endpoint on triangle plane\n") end  
    if (sign(s0) * sign(s1) >= 0) return false end

    # Volumes with each triangle edge; all must share sign for interior hit
    t0 = signed_volume6(p0, p1, v0, v1)
    t1 = signed_volume6(p0, p1, v1, v2)
    t2 = signed_volume6(p0, p1, v2, v0)

    pos = (t0 > 0 && t1 > 0 && t2 > 0)
    neg = (t0 < 0 && t1 < 0 && t2 < 0)
    return (pos || neg)
end

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
function t_index(a::Int, b::Int, c::Int)
    mn = min(a,b,c)
    mx = max(a,b,c)
    md = a-mn+c-mx+b
    return triangle_index(mn, md, mx)  # +1 for 1-based indexing
end

@inline function edge_index(a::Int, b::Int)
    @boundscheck (1 <= a < b <= 16) || throw(ArgumentError("edge endpoints must satisfy 1 <= a < b <= 16"))
    return a + (b-1)*(b-2) ÷ 2 # 1 based
end
e_index(a::Int, b::Int) = @inbounds edge_index(minmax(a,b)...) # 1 based
@propagate_inbounds singleton(a::Int, b::Int) = BitSet128.singleton(edge_index(a,b))

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

# ==== Small isbits structs ====

struct VertexCorner
    next::Int8
    opp::Int8
    tri::Int16
end

struct VertexState
    avail_tri::Int8
    loc::Int8
    min_tri::Int16
end

# ==== Stack implementation ====

mutable struct SmallStack{T}
    st::Vector{T}
    sp::Int
end

@inline function SmallStack(default::T, nstack::Int) where {T}
    st = Vector{T}(undef, nstack)
    fill!(st, default)
    return SmallStack{T}(st, 0)
end

@inline isempty(s::SmallStack) = (s.sp == 0)

@inline function Base.push!(s::SmallStack{T}, x::T) where {T}
    @boundscheck s.sp < length(s.st) || throw(BoundsError())
    sp = s.sp + 1
    s.sp = sp
    st = s.st
    @inbounds st[sp] = x
    return s
end

@inline function Base.pop!(s::SmallStack{T}) where {T}
    @boundscheck s.sp > 0 || throw(ArgumentError("empty stack"))
    st = s.st
    sp = s.sp
    @inbounds x = st[sp]
    s.sp = sp - 1
    return x
end

@inline peek(s::SmallStack{T}) where {T} = (@inbounds s.st[s.sp])

# ==== Factory ====

@inline function make_list_stack(default_obj::T, nlist::Int, nstack::Int) where {T}
    list = Vector{T}(undef, nlist)
    fill!(list, default_obj)
    stack = SmallStack(default_obj, nstack)
    return list, stack
end

# ==== Example usage ====

# VertexCorner storage + stack
vc, cs = make_list_stack(VertexCorner(0,0,0), 16, 44)

# VertexState storage + stack
vstate, vst = make_list_stack(
    VertexState(edge_index(1,2), 0, triangle_index(15,15,16)),
    16, 16
)

# ==== Quick sanity checks ====

@assert isbitstype(VertexCorner)
@assert isbitstype(VertexState)

# Should be zero allocations
@assert @allocated(begin
    push!(cs, vc[1])
    pop!(cs)
end) == 0



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
    results = Surface[]
    for t in min_tri:max_tri
        tri_map, ETConf, curr_edgeset, eset, te_conf, ti_map = initialize_with_maxtriangle(max_tri, conflictcount, edgeset, triangle_map, ETConflict, points)
    
    
    return results
end
end


end # module


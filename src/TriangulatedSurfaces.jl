module TriangulatedSurfaces

using Random

export Point3D, Triangle, Surface, 
       tetrahedron_with_origin, add_random_interior_points!,
       has_colinear_triple, has_coplanar_quad, edge_index, triangle_index, conflict,
       enumerate_triangulated_surfaces, precompute_maps

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
@inline crossp(u::NTuple{3,<:Integer}, v::NTuple{3,<:Integer}) = (
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
    has_colinear_triple(points; atol=1e-9)

Return `true` if any triple of points is colinear (within `atol`).
"""
function has_colinear_triple(points::Vector{Point3D}; atol::Int=0)
    n = length(points)
    n < 3 && return false
    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        u = sub(points[j], points[i])
        v = sub(points[k], points[i])
        if norm2(crossp(u, v)) <= atol
            return true
        end
    end
    return false
end

"""
    has_coplanar_quad(points; atol=1e-9)

Return `true` if any set of four points is coplanar (within `atol`).
"""
function has_coplanar_quad(points::Vector{Point3D})
    n = length(points)
    n < 4 && return false
    for i in 1:(n - 3), j in (i + 1):(n - 2), k in (j + 1):(n - 1), l in (k + 1):n
        u = sub(points[j], points[i])
        v = sub(points[k], points[i])
        w = sub(points[l], points[i])
        if dotp(crossp(u, v), w) == 0
            return true
        end
    end
    return false
end

"""
    point_in_tetra(p, v1, v2, v3, v4; atol=1e-9)

Return `true` if `p` lies strictly inside the tetrahedron defined by
vertices `v1`..`v4`. Points on the boundary are rejected.
"""
function point_in_tetra(p::Point3D, v1::Point3D, v2::Point3D, v3::Point3D, v4::Point3D; atol::Int=0)
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
    edge_index(a, b)

Map an edge with endpoints `a < b` to an integer: `a + binomial(b-1, 2)`.
"""
edge_index(a::Int, b::Int) = a + (b-1)*(b-2) ÷ 2


"""
    triangle_index(a, b, c)

Map a triangle `a < b < c` to an integer: `a + binomial(b-1, 2) + binomial(c-1, 3)`.
"""
triangle_index(a::Int, b::Int, c::Int) = a + (b-1)*(b-2) ÷ 2 + (c-1)*(c-2)*(c-3) ÷ 6



"""
    precompute_maps(points::Vector{Point3D})

Precompute mappings:
1. edge_map: map from edge numbers 1:n choose 2 to pairs of vertex indices (a,b) with a<b
1. triangle_map: map from triangle numbers 1:n choose 3 to triples of vertex indices (c,d,e) with c<d<e
1. ETConflict: Boolean matrix where ETConflict[edge_index(a,b), triangle_index(c,d,e)] 
    is true iff edge (a,b) conflicts with triangle (c,d,e)
2. triangle_edges: (n choose 3)×3 matrix where row triangle_index(a,b,c) contains
    [edge_index(a,b), edge_index(b,c), edge_index(a,c)]
3. edge_triangles: (n choose 2)×(n-2) matrix where row edge_index(a,b) contains
    triangle_index values for all triangles containing that edge
"""
function precompute_maps(points::Vector{Point3D})
    n = length(points)
    
    # Maximum possible indices
    max_edge_idx = edge_index(n - 1, n)
    max_tri_idx = triangle_index(n - 2, n - 1, n)
    
    # Create reverse maps: index -> edge and index -> triangle
    edge_map = Vector{Tuple{Int,Int}}(undef, max_edge_idx)
    triangle_map = Vector{NTuple{3,Int}}(undef, max_tri_idx)

    for a in 1:(n-1), b in (a+1):n
        edge_map[edge_index(a, b)] = (a, b)
    end

    for c in 1:(n-2), d in (c+1):(n-1), e in (d+1):n
        triangle_map[triangle_index(c, d, e)] = (c, d, e)
    end
    
    # Initialize Boolean matrix for conflicts
    ETConflict = falses(max_edge_idx, max_tri_idx)
    
    # Initialize triangle_edges: each row is 3 edge indices
    triangle_edges = zeros(Int, max_tri_idx, 3)
    
    # Populate triangle_edges first (needed for other computations)
    for tri_idx in 1:max_tri_idx
        c, d, e = triangle_map[tri_idx]
        triangle_edges[tri_idx, 1] = edge_index(c, d)
        triangle_edges[tri_idx, 2] = edge_index(d, e)
        triangle_edges[tri_idx, 3] = edge_index(c, e)
    end
    
    # Initialize edge_triangles: list of triangles containing each edge
    edge_triangles = [Int[] for _ in 1:max_edge_idx]
    for tri_idx in 1:max_tri_idx
        for j in 1:3
            edge_idx = triangle_edges[tri_idx, j]
            push!(edge_triangles[edge_idx], tri_idx)
        end
    end
    
    # Populate conflict matrix: loop over edge and triangle numbers
    for edge_idx in 1:max_edge_idx
        a, b = edge_map[edge_idx]
        for tri_idx in 1:max_tri_idx
            c, d, e = triangle_map[tri_idx]
            ETConflict[edge_idx, tri_idx] = conflict((a, b), (c, d, e), points)
        end
    end
    
    # Build conflict lists from ETConflict matrix
    # edge_conflicts[e] = list of triangle indices that edge e conflicts with
    edge_conflicts = [Int[] for _ in 1:max_edge_idx]
    # triangle_conflicts[t] = list of edge indices that triangle t conflicts with
    triangle_conflicts = [Int[] for _ in 1:max_tri_idx]
    
    for edge_idx in 1:max_edge_idx
        for tri_idx in 1:max_tri_idx
            if ETConflict[edge_idx, tri_idx]
                push!(edge_conflicts[edge_idx], tri_idx)
                push!(triangle_conflicts[tri_idx], edge_idx)
            end
        end
    end
    
    return edge_map, triangle_map, triangle_edges, edge_triangles, edge_conflicts, triangle_conflicts
end





# --- Internal helpers ---


# Segment–triangle intersection using only integer signed volumes.
# Assumes non-degenerate input (no coplanar quadruples).
@inline signed_volume6(a::Point3D, b::Point3D, c::Point3D, d::Point3D) = dotp(crossp(sub(b, a), sub(c, a)), sub(d, a))

function segment_intersects_triangle(p0::Point3D, p1::Point3D, v0::Point3D, v1::Point3D, v2::Point3D)
    # Endpoints on opposite sides of the triangle plane
    s0 = signed_volume6(v0, v1, v2, p0)
    s1 = signed_volume6(v0, v1, v2, p1)
    if (s0 == 0 || s1 == 0) print("Warning: segment endpoint on triangle plane\n") end  
    (s0 == 0 || s1 == 0 || s0 * s1 > 0) && return false

    # Volumes with each triangle edge; all must share sign for interior hit
    t0 = signed_volume6(p0, p1, v0, v1)
    t1 = signed_volume6(p0, p1, v1, v2)
    t2 = signed_volume6(p0, p1, v2, v0)

    pos = (t0 > 0 && t1 > 0 && t2 > 0)
    neg = (t0 < 0 && t1 < 0 && t2 < 0)
    return (pos || neg)
end

# --- Triangulated Surface Enumeration via Backtracking ---

"""
    enumerate_triangulated_surfaces(points; max_results=1000)

Enumerate all valid triangulated surfaces on `points` using backtracking.
Returns a vector of `Surface` instances satisfying:
1. Manifold: every edge in 0 or 2 triangles
2. No intersection: edges don't intersect triangles except at shared vertices
3. Connected: path exists between any two vertices
4. Link: triangles around each vertex form a single loop

Limited to `max_results` surfaces to prevent exhaustive search.
"""
function enumerate_triangulated_surfaces(points::Vector{Point3D}; max_results::Int=1000)
    n = length(points)
    n < 4 && throw(ArgumentError("Need at least 4 points"))
    
    # Precompute conflict matrix and mappings
    edge_map, triangle_map, triangle_edges, edge_triangles, edge_conflicts, triangle_conflicts = precompute_maps(points)
    
    # Generate all possible triangle indices
    all_tri_indices = 1:triangle_index(n - 2, n - 1, n)
    
    results = Surface[]
    state = BacktrackState(n, points, edge_map, triangle_map, triangle_edges)
    
    backtrack!(state, all_tri_indices, edge_conflicts, triangle_conflicts, 1, results, max_results)
    return results
end

mutable struct BacktrackState
    n::Int                              # number of points
    points::Vector{Point3D}             # the point set
    triangle_indices::Vector{Int}       # currently selected triangle indices
    edge_count::Vector{Int}             # edge_index -> count (0, 1, or 2)
    edge_indices::Set{Int}              # all edge indices used
    edge_map::Vector{Tuple{Int,Int}}    # index -> (a, b) edge
    triangle_map::Vector{NTuple{3,Int}} # index -> (c, d, e) triangle
    triangle_edges::Matrix{Int}         # tri_idx -> [edge1, edge2, edge3]
end

function BacktrackState(n::Int, points::Vector{Point3D}, edge_map::Vector{Tuple{Int,Int}}, 
                       triangle_map::Vector{NTuple{3,Int}}, triangle_edges::Matrix{Int})
    max_edge_idx = edge_index(n - 1, n)
    BacktrackState(n, points, Int[], zeros(Int, max_edge_idx), Set{Int}(), edge_map, triangle_map, triangle_edges)
end

function add_triangle!(state::BacktrackState, tri_idx::Int)
    push!(state.triangle_indices, tri_idx)
    for j in 1:3
        edge_idx = state.triangle_edges[tri_idx, j]
        state.edge_count[edge_idx] += 1
        push!(state.edge_indices, edge_idx)
    end
end

function remove_triangle!(state::BacktrackState, tri_idx::Int)
    pop!(state.triangle_indices)
    for j in 1:3
        edge_idx = state.triangle_edges[tri_idx, j]
        state.edge_count[edge_idx] -= 1
        if state.edge_count[edge_idx] == 0
            delete!(state.edge_indices, edge_idx)
        end
    end
end

function backtrack!(state::BacktrackState, all_tri_indices, edge_conflicts::Vector{Vector{Int}},
                   triangle_conflicts::Vector{Vector{Int}}, idx::Int, results::Vector{Surface}, max_results::Int)
    length(results) >= max_results && return
    
    # Try completing surface with current triangles
    if is_valid_surface(state)
        triangles = [Triangle(state.triangle_map[ti]...) for ti in state.triangle_indices]
        surface = Surface(state.points, triangles)
        push!(results, surface)
        return
    end
    
    # Try adding more triangles
    idx > length(all_tri_indices) && return
    
    for i in idx:length(all_tri_indices)
        tri_idx = all_tri_indices[i]
        
        # Quick pruning checks before adding
        if can_add_triangle(state, tri_idx, edge_conflicts, triangle_conflicts)
            add_triangle!(state, tri_idx)
            backtrack!(state, all_tri_indices, edge_conflicts, triangle_conflicts, i + 1, results, max_results)
            remove_triangle!(state, tri_idx)
        end
    end
end

function can_add_triangle(state::BacktrackState, tri_idx::Int, 
                         edge_conflicts::Vector{Vector{Int}}, triangle_conflicts::Vector{Vector{Int}})::Bool
    # Check manifold constraint: no edge would exceed count 2
    for j in 1:3
        edge_idx = state.triangle_edges[tri_idx, j]
        if state.edge_count[edge_idx] >= 2
            return false
        end
    end
    
    # Check intersection using precomputed conflict lists:
    # new triangle's edges don't conflict with existing triangles
    for j in 1:3
        edge_idx = state.triangle_edges[tri_idx, j]
        for existing_tri_idx in edge_conflicts[edge_idx]
            if existing_tri_idx in state.triangle_indices
                return false
            end
        end
    end
    
    # Check existing edges don't conflict with new triangle
    for existing_edge_idx in triangle_conflicts[tri_idx]
        if existing_edge_idx in state.edge_indices
            return false
        end
    end
    
    return true
end

function is_valid_surface(state::BacktrackState)::Bool
    isempty(state.triangle_indices) && return false
    
    # 1. Manifold: all edges in exactly 2 triangles
    for edge_idx in state.edge_indices
        if state.edge_count[edge_idx] != 2
            return false
        end
    end
    
    # 2. No intersection already checked during construction
    
    # 3. Connected: all vertices reachable
    is_connected(state) || return false
    
    # 4. Link: triangles around each vertex form a loop
    check_vertex_links(state) || return false
    
    return true
end

function is_connected(state::BacktrackState)::Bool
    isempty(state.edge_indices) && return false
    
    # Find vertices that appear in the triangulation
    vertices = Set{Int}()
    for tri_idx in state.triangle_indices
        c, d, e = state.triangle_map[tri_idx]
        union!(vertices, [c, d, e])
    end
    isempty(vertices) && return false
    
    # BFS from first vertex
    start = minimum(vertices)
    visited = Set{Int}([start])
    queue = [start]
    
    while !isempty(queue)
        v = popfirst!(queue)
        for edge_idx in state.edge_indices
            u, w = state.edge_map[edge_idx]
            neighbor = (u == v) ? w : (w == v) ? u : 0
            if neighbor > 0 && !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    return length(visited) == length(vertices)
end

function check_vertex_links(state::BacktrackState)::Bool
    # Build adjacency: vertex -> list of triangle indices containing it
    vertex_triangles = Dict{Int, Vector{Int}}()
    
    for tri_idx in state.triangle_indices
        c, d, e = state.triangle_map[tri_idx]
        for v in [c, d, e]
            if !haskey(vertex_triangles, v)
                vertex_triangles[v] = Int[]
            end
            push!(vertex_triangles[v], tri_idx)
        end
    end
    
    # Check each vertex has a valid link (forms a single loop)
    for (v, tri_indices) in vertex_triangles
        if !is_valid_link(v, tri_indices, state)
            return false
        end
    end
    
    return true
end

function is_valid_link(vertex::Int, tri_indices::Vector{Int}, state::BacktrackState)::Bool
    isempty(tri_indices) && return false
    
    # Extract edges around this vertex
    edges = Tuple{Int,Int}[]
    for tri_idx in tri_indices
        c, d, e = state.triangle_map[tri_idx]
        if c == vertex
            push!(edges, (min(d, e), max(d, e)))
        elseif d == vertex
            push!(edges, (min(c, e), max(c, e)))
        else  # e == vertex
            push!(edges, (min(c, d), max(c, d)))
        end
    end
    
    # Build adjacency for the link
    link_adj = Dict{Int, Vector{Int}}()
    for (u, v) in edges
        if !haskey(link_adj, u)
            link_adj[u] = Int[]
        end
        if !haskey(link_adj, v)
            link_adj[v] = Int[]
        end
        push!(link_adj[u], v)
        push!(link_adj[v], u)
    end
    
    # Check each vertex in link has exactly 2 neighbors (forms a cycle)
    for (_, neighbors) in link_adj
        length(neighbors) != 2 && return false
    end
    
    # Check it's a single cycle (connected)
    isempty(link_adj) && return false
    start = first(keys(link_adj))
    visited = Set{Int}([start])
    current = start
    prev = -1
    
    for _ in 1:length(link_adj)
        neighbors = link_adj[current]
        next = (neighbors[1] == prev) ? neighbors[2] : neighbors[1]
        next == start && break
        next in visited && return false  # cycle too short
        push!(visited, next)
        prev = current
        current = next
    end
    
    return length(visited) == length(link_adj)
end

end # module


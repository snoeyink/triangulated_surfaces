"""
    Point3D(x, y, z)

Minimal 3D point type.
"""
struct Point3D
    x::Int
    y::Int
    z::Int
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

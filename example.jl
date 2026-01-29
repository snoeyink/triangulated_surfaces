include(joinpath(@__DIR__, "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

using Random

println("Step 1: Loading module... done")

# Start with tetrahedron vertices plus the origin.
println("Step 2: Creating tetrahedron...")
tetra = tetrahedron_with_origin()


# Add up to 6 random interior integer points.
println("Step 3: Sampling random interior points...")
rng = MersenneTwister(25122025)
all_points = add_random_interior_points!(tetra, 9; rng, max_attempts=10000)

println("\nStep 4: Basic checks")
println("Total points: ", length(all_points))
println("Colinear triple? ", has_colinear_triple(all_points))
println("Coplanar quad? ", has_coplanar_quad(all_points))

# Example edge/triangle ids and conflict test.
println("\nStep 5: Index and conflict test")
edge = (1, 5)
tri = (2, 3, 6)
println("edge index: ", edge_index(edge...))
println("triangle index: ", triangle_index(tri...))
println("conflict(edge, tri): ", conflict(edge, tri, all_points))

# Precompute conflict maps
println("\nStep 6: Precomputing conflict maps...")
edge_map, tri_map, triangle_edges, edge_triangles, edge_conflicts, triangle_conflicts = precompute_maps(all_points)

println("Edge map length: ", length(edge_map))
println("Triangle map length: ", length(tri_map))
println("Triangle edges size: ", size(triangle_edges))
println("Edge triangles length: ", length(edge_triangles))
println("Edge conflicts length: ", length(edge_conflicts))
println("Triangle conflicts length: ", length(triangle_conflicts))
println("Example: edge_map[$(edge_index(edge...))] = ", edge_map[edge_index(edge...)])
println("Example: tri_map[$(triangle_index(tri...))] = ", tri_map[triangle_index(tri...)])

println("\nStep 7: Conflict length histograms")
# Histogram of edge_conflicts lengths
edge_len_counts = Dict{Int,Int}()
for lst in edge_conflicts
	k = length(lst)
	edge_len_counts[k] = get(edge_len_counts, k, 0) + 1
end
println("Edge conflict lengths (len => count):")
for k in sort(collect(keys(edge_len_counts)))
	println("  ", k, " => ", edge_len_counts[k])
end

# Histogram of triangle_conflicts lengths
tri_len_counts = Dict{Int,Int}()
for lst in triangle_conflicts
	k = length(lst)
	tri_len_counts[k] = get(tri_len_counts, k, 0) + 1
end
println("Triangle conflict lengths (len => count):")
for k in sort(collect(keys(tri_len_counts)))
	println("  ", k, " => ", tri_len_counts[k])
end

println("\nDone!")

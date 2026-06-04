using GLMakie
GLMakie.activate!()

const ueindex = TriangulatedSurfaces.ueindex

# Plot state structure
mutable struct TriPlot
    fig::Figure
    ax::Axis3
    vertices::Vector{GLMakie.Point3f}
    triangles::Vector{NamedTuple{(:indices, :mesh, :edges, :labels), 
                                  Tuple{Tuple{Int,Int,Int}, Mesh, Vector{Lines}, Vector{GLMakie.Text}}}}
end

function create_plot(points::Vector{Point3D})
    @assert length(points) <= 16 "Maximum 16 points allowed"

    vertices = GLMakie.Point3f[(p.x, p.y, p.z) for p in points]
    
    fig = Figure(size=(1200, 1000))
    
    # Compute bounding box
    xs = [p[1] for p in vertices]
    ys = [p[2] for p in vertices]
    zs = [p[3] for p in vertices]
    
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    zmin, zmax = extrema(zs)
    
    # Add 10% margin
    margin = 0.1
    xrange = xmax - xmin
    yrange = ymax - ymin
    zrange = zmax - zmin
    
    ax = Axis3(fig[1, 1],
               xlabel="X", ylabel="Y", zlabel="Z",
               title="Triangulated Surface",
               aspect=:data,
               limits=(xmin - margin*xrange, xmax + margin*xrange,
                      ymin - margin*yrange, ymax + margin*yrange,
                      zmin - margin*zrange, zmax + margin*zrange))
    
    # Plot points as spheres
    scatter!(ax, xs, ys, zs, 
             color=:blue, 
             markersize=15)
    
    # Add point labels
    for (i, p) in enumerate(vertices)
        text!(ax, string(i),
              position=(p[1], p[2], p[3]),
              fontsize=20,
              color=:black,
              align=(:center, :center),
              offset=(0, 0, 10))
    end
    
    display(fig)
    
    return TriPlot(fig, ax, vertices, [])
end

function plot_tri(pl::TriPlot, a::Integer, b::Integer, c::Integer)
    @assert 1 <= a <= length(pl.vertices) "Point index $a out of range"
    @assert 1 <= b <= length(pl.vertices) "Point index $b out of range"
    @assert 1 <= c <= length(pl.vertices) "Point index $c out of range"
    @assert allunique([a, b, c]) "Triangle vertices must be distinct"
    
    # Check for duplicate triangle
    for tri in pl.triangles
        if Set(tri.indices) == Set((a, b, c))
            @warn "Triangle ($a, $b, $c) already exists!"
            return
        end
    end
    
    # Change previous triangle from red to green if it exists
    if !isempty(pl.triangles)
        last_tri = pl.triangles[end]
        last_tri.mesh.color = RGBAf(0, 1, 0, 0.3)
    end
    
    # Get points
    pa = pl.vertices[a]
    pb = pl.vertices[b]
    pc = pl.vertices[c]
    
    # Create triangle mesh (new one is red)
    vertices = GLMakie.Point3f[pa, pb, pc]
    faces = GLMakie.GLTriangleFace[(1, 2, 3)]
    tri_mesh = mesh!(pl.ax, vertices, faces,
                     color=RGBAf(1, 0, 0, 0.3),  # Red, semi-transparent
                     transparency=true)
    
    # Draw edges and labels
    edges_plot = Lines[]
    labels_plot = GLMakie.Text[]
    
    for (i1, i2) in [(a, b), (b, c), (c, a)]
        p1 = pl.vertices[i1]
        p2 = pl.vertices[i2]
        
        # Draw edge
        edge_line = lines!(pl.ax, 
                          [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]],
                          color=:black,
                          linewidth=2)
        push!(edges_plot, edge_line)
        
        # Edge label at midpoint
        mid_x = (p1[1] + p2[1]) / 2
        mid_y = (p1[2] + p2[2]) / 2
        mid_z = (p1[3] + p2[3]) / 2
        
        edge_label = text!(pl.ax, string(ueindex(i1, i2)),
                          position=(mid_x, mid_y, mid_z),
                          fontsize=14,
                          color=:red,
                          align=(:center, :center))
        push!(labels_plot, edge_label)
    end
    
    # Store triangle data
    push!(pl.triangles, (indices=(a, b, c), 
                         mesh=tri_mesh, 
                         edges=edges_plot, 
                         labels=labels_plot))
    
    return nothing
end

function plot_pop!(pl::TriPlot)
    if isempty(pl.triangles)
        @warn "No triangles to remove"
        return
    end
    
    # Remove last triangle
    tri = pop!(pl.triangles)
    
    # Delete mesh
    delete!(pl.ax, tri.mesh)
    
    # Delete edges
    for edge in tri.edges
        delete!(pl.ax, edge)
    end
    
    # Delete labels
    for label in tri.labels
        delete!(pl.ax, label)
    end
    
    # If there are remaining triangles, make the new last one red
    if !isempty(pl.triangles)
        pl.triangles[end].mesh.color = RGBAf(1, 0, 0, 0.3)
    end
    
    return nothing
end
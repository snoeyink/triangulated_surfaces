using Test

include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
using .TriangulatedSurfaces

include("bitset.jl")

@testset "TriangulatedSurfaces basic geometry" begin
    points = tetrahedron_with_origin(scale=4)

    @test points == [
        Point3D(4, 4, 4),
        Point3D(4, -4, -4),
        Point3D(-4, 4, -4),
        Point3D(-4, -4, 4),
        Point3D(0, 0, 0),
    ]

    @test TriangulatedSurfaces.point_in_tetra(Point3D(-2, -2, -2), points[1:4]...) == false
    @test TriangulatedSurfaces.point_in_tetra(Point3D(-1, -1, -1), points[1:4]...) == true
    @test TriangulatedSurfaces.point_in_tetra(Point3D(2, 1, 0), points[1:4]...) == true

    push!(points, Point3D(-1, -1, -1))
    @test has_coplanar_quad(points) == true

    pop!(points)
    push!(points, Point3D(1, 0, 0))
    @test has_coplanar_quad(points) == true

    pop!(points)
    push!(points, Point3D(2, 1, 0))
    @test has_coplanar_quad(points) == false
end

@testset "TriangulatedSurfaces conflicts" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,et,es = TriangulatedSurfaces.precompute_conflicts(points)
    @test Tuple.(findall(et)) == [(TriangulatedSurfaces.edge_index(4,6),19), (TriangulatedSurfaces.edge_index(1,5),20)]
    @test tm[19] == (1, 2, 5)
    @test tm[20] == (3, 4, 6)
    @test es[19] == TriangulatedSurfaces.singleton(4,6)
    @test es[20] == TriangulatedSurfaces.singleton(1,5)
    @test isdisjoint(es[19], es[20]) == true
    @test isdisjoint(es[19], TriangulatedSurfaces.singleton(4,6)|TriangulatedSurfaces.singleton(1,5)) == false
end


@testset "TriangulatedSurfaces boundary loop" begin
    vc, cs = TriangulatedSurfaces.make_list_stack(vertex_corner(0,0,0), 16, 44)
    vstate, vst = TriangulatedSurfaces.make_list_stack(vertex_state(0, triangle_index(15,15,16)), 16, 16)
    
    nbdry = 6
    vc[1] = TriangulatedSurfaces.vertex_corner(3, 5, TriangulatedSurfaces.triangle_index(1,3,5))
    vc[3] = TriangulatedSurfaces.vertex_corner(6, 2, TriangulatedSurfaces.triangle_index(2,3,6))
    vc[6] = TriangulatedSurfaces.vertex_corner(4, 7, TriangulatedSurfaces.triangle_index(4,6,7))
    vc[4] = TriangulatedSurfaces.vertex_corner(7, 6, TriangulatedSurfaces.triangle_index(4,6,7))
    vc[7] = TriangulatedSurfaces.vertex_corner(5, 6, TriangulatedSurfaces.triangle_index(5,6,7))
    vc[5] = TriangulatedSurfaces.vertex_corner(1, 3, TriangulatedSurfaces.triangle_index(1,3,5))
    vc[2] = TriangulatedSurfaces.vertex_corner(6, 5, TriangulatedSurfaces.triangle_index(2,5,6)) # now interior
    TriangulatedSurfaces.push!(cs, vc[2])
    reme = 5
    
    vc[11] = TriangulatedSurfaces.vertex_corner(12, 13, TriangulatedSurfaces.triangle_index(11,12,13))
    vc[12] = TriangulatedSurfaces.vertex_corner(13, 11, TriangulatedSurfaces.triangle_index(11,12,13))
    vc[13] = TriangulatedSurfaces.vertex_corner(11, 12, TriangulatedSurfaces.triangle_index(11,12,13))

    e1(a) = @inbounds begin b = vc[a].next; c = vc[a].opp; return vc[b].next == c end
    e2(a) = @inbounds begin c = vc[a].opp; return vc[c].next == a end
    @test e1.([1,3,6,4,7,5,11,12,13,2]) == [false, false, true, false, false, true, true, true, true, false]
    @test e2.([1,3,6,4,7,5,11,12,13,2]) == [true, false, false, true, false, false, true, true, true, false]
  
    function removable(vc, a, nbdry)   
        bv = falses(16)
        for i = 1:nbdry
            bv[a] = true
            a = vc[a].next
        end
        mn = TriangulatedSurfaces.triangle_index(14,15,16)+1 # larger than any triangle index
        mi = a
        for i = 1:nbdry
            if vc[a].tri < mn && (!bv[vc[a].opp] || e1(a))
               mn = vc[a].tri
               mi = a
            end
            a = vc[a].next
        end
        return mn, mi
    end
    @test removable(vc, 1, nbdry)[2] == reme
    @test removable(vc, 12, 3)[2] == 12
end

@testset "TriangulatedSurfaces vertex state and boundary vertices" begin
    vc, cs = TriangulatedSurfaces.make_list_stack(vertex_corner(0,0,0), 16, 44)
    vstate, vst = TriangulatedSurfaces.make_list_stack(vertex_state(0, triangle_index(15,15,16)), 16, 16)
end

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


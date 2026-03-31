
@testset "TriangulatedSurfaces conflicts" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,es = TriangulatedSurfaces.precompute_conflicts(points)
    @test tm[19] == (1, 2, 5)
    @test tm[20] == (3, 4, 6)
    @test es[19].conf == TriangulatedSurfaces.singleton(4,6)
    @test es[20].conf == TriangulatedSurfaces.singleton(1,5)
    @test isdisjoint(es[19].conf, es[20].conf) == true
    @test !isdisjoint(es[19].conf, TriangulatedSurfaces.singleton(4,6)|TriangulatedSurfaces.singleton(1,5))
end


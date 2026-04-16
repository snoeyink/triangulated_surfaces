const singleton = TriangulatedSurfaces.singleton
const precompute_conflicts = TriangulatedSurfaces.precompute_conflicts
const build_tri_table = TriangulatedSurfaces.build_tri_table
@testset "TriangulatedSurfaces conflicts" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,es = precompute_conflicts(points)
    print(tm)

    @test tm[19] == (1, 2, 5)
    @test tm[20] == (3, 4, 6)
    @test es[19].conf == singleton(4,6)
    @test es[20].conf == singleton(1,5)
    @test isdisjoint(es[19].conf, es[20].conf) 
    @test !isdisjoint(es[19].conf, singleton(4,6)|singleton(1,5))
end

@testset "build_tri_table" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,es = precompute_conflicts(points)
    tmax, tmap, esets, tri_table = build_tri_table(UInt8(6), UInt16(length(es)), tm, es)
    #println(tmap)
    @test tmax == 16 # 16 triangles survive 
    @test tri_table[3, 4, 6] == tmax
    @test tri_table[1, 2, 5] == tmax+1  # conflicts with tmax, so sentinel
    @test tri_table[1, 2, 3] == 1     # no conflict, so renumbered index is just original index
    @test tri_table[1, 4, 5] == tmax+1  # conflict with edge 1,5, so tmax
    @test tri_table[1, 3, 4] == 3     # no conflict, so renumbered index is just original index
    @test tri_table[1, 2, 6] == 8    # no conflicts with tmax

end

@testset "build_tri_table2" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,es = precompute_conflicts(points)
    tmax, tmap, esets, tri_table = build_tri_table(6, length(es)-1, tm, es)
    #println(tmap)
    @test tmax == 16 
    @test tri_table[3, 4, 6] == tmax+1 # conflicts with tmax, so sentinel
    @test tri_table[1, 2, 5] == tmax # is tmax. 
    @test tri_table[1, 2, 3] == 1     # no conflict, so renumbered index is just original index
    @test tri_table[1, 4, 5] == 14  # conflict with edge 1,5, so tmax
    @test tri_table[1, 3, 4] == 3     # no conflict, so renumbered index is just original index
    @test tri_table[4, 5, 6] == tmax+1    # no conflicts with tmax

end
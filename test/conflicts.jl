const singleton = TriangulatedSurfaces.singleton
const precompute_conflicts = TriangulatedSurfaces.precompute_conflicts
const triangle_index = TriangulatedSurfaces.triangle_index
const edge_index = TriangulatedSurfaces.edge_index
const BitSet128 = TriangulatedSurfaces.BitSet128

@testset "TriangulatedSurfaces indices and conflicts" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm, edge_map, econfl = precompute_conflicts(points)

    @test length(edge_map) == edge_index(length(points), length(points) - 1)

    for a in 1:length(points), b in 1:length(points)
        a == b && continue
        @test edge_map[edge_index(a, b)] == (a, b)
    end

    for a in 1:length(points), b in (a + 1):length(points)
        @test edge_index(b, a) == edge_index(a, b) + 128
    end

    for a in 1:length(points), b in 1:length(points), c in 1:length(points)
        (a == b || b == c || a == c) && continue
        tri = tm[triangle_index(a, b, c)]
        @test tri[1] != 0
        @test triangle_index(tri...) == triangle_index(a, b, c)
    end

    t125 = triangle_index(1, 2, 5)
    t346 = triangle_index(3, 4, 6)
    @test tm[t125][1] != 0
    @test tm[t346][1] != 0
    @test triangle_index(tm[t125]...) == t125
    @test triangle_index(tm[t346]...) == t346


    @test tm[t125] == (1, 2, 5)
    @test tm[t346] == (3, 4, 6)

    expected_125 = singleton(BitSet128, edge_index(4, 6))
    expected_346 = singleton(BitSet128, edge_index(1, 5))
    expected_union = expected_125 | expected_346

    @test econfl[t125] == expected_125
    @test econfl[t346] == expected_346
    @test isdisjoint(econfl[t125], econfl[t346])
    @test !isdisjoint(econfl[t125], expected_union)
    @test !isdisjoint(econfl[t346], expected_union)
end

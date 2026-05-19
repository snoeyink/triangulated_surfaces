using Test

if !isdefined(@__MODULE__, :TriangulatedSurfaces)
    include(joinpath(@__DIR__, "..", "src", "TriangulatedSurfaces.jl"))
    using .TriangulatedSurfaces
end

const PUF = TriangulatedSurfaces.PackedUnionFind

@testset "PackedUnionFind" begin
    @testset "Initial roots" begin
        uf = PUF.PackedUF()
        for i in 1:16
            @test PUF.find_root(uf, i) == i
        end
    end

    @testset "Basic unions" begin
        uf0 = PUF.PackedUF()

        uf1, ok1 = PUF.union_sets(uf0, 1, 2)
        @test ok1
        @test PUF.find_root(uf1, 1) == 1
        @test PUF.find_root(uf1, 2) == 1

        uf2, ok2 = PUF.union_sets(uf1, 2, 3)
        @test ok2
        @test PUF.find_root(uf2, 1) == 1
        @test PUF.find_root(uf2, 2) == 1
        @test PUF.find_root(uf2, 3) == 1

        @test PUF.find_root(uf2, 4) == 4
    end

    @testset "Cycle detection and immutability" begin
        uf0 = PUF.PackedUF()
        uf1, ok1 = PUF.union_sets(uf0, 4, 5)
        @test ok1

        uf2, ok2 = PUF.union_sets(uf1, 5, 6)
        @test ok2

        uf3, ok3 = PUF.union_sets(uf2, 4, 6)
        @test !ok3
        @test uf3 === uf2

        @test PUF.find_root(uf2, 4) == PUF.find_root(uf2, 5)
        @test PUF.find_root(uf2, 5) == PUF.find_root(uf2, 6)
    end

    @testset "Merge existing components" begin
        uf0 = PUF.PackedUF()

        uf1, ok1 = PUF.union_sets(uf0, 1, 2)
        uf2, ok2 = PUF.union_sets(uf1, 3, 4)
        uf3, ok3 = PUF.union_sets(uf2, 2, 4)

        @test ok1
        @test ok2
        @test ok3

        r1 = PUF.find_root(uf3, 1)
        r2 = PUF.find_root(uf3, 2)
        r3 = PUF.find_root(uf3, 3)
        r4 = PUF.find_root(uf3, 4)

        @test r1 == r2 == r3 == r4
    end
end

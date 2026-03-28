@testset "BitSet128" begin
    @test isbitstype(BitSet128)

    empty = BitSet128()
    @test empty.bits == zero(UInt128)
    @test iszero(empty)
    @test length(empty) == 0

    low = TriangulatedSurfaces.singleton(1)
    mid = TriangulatedSurfaces.singleton(64)
    high = TriangulatedSurfaces.singleton(128)
    pair = TriangulatedSurfaces.push(low, 128)
    @test (1 in pair) && (128 in pair) 

    @test 1 in low
    @test !(2 in low)
    @test 64 in mid
    @test 128 in high
    @test length(TriangulatedSurfaces.push(low, 1)) == 1

    @test (low & high).bits == zero(UInt128)
    @test (low | high).bits == pair.bits
    @test (pair ⊻ low).bits == high.bits 


    complement = ~low
    @test !(1 in complement)
    @test 2 in complement
    @test 128 in complement
    @test length(complement) == 127

    @test issubset(low, pair)
    @test issubset(high, pair)
    @test !issubset(pair, low)

    a = BitSet128(UInt128(0b1010))
    b = BitSet128(UInt128(0b1100))
    @test (a & b).bits == UInt128(0b1000)
    @test (a | b).bits == UInt128(0b1110)
    @test ⊻(a, b).bits == UInt128(0b0110) 

    @test_throws ArgumentError (0 in empty)
    @test_throws ArgumentError (129 in empty)
    @test_throws ArgumentError TriangulatedSurfaces.push(empty, 0)
    @test_throws ArgumentError TriangulatedSurfaces.push(empty, 129)
    @test_throws ArgumentError Base.delete(empty, 0)
    @test_throws ArgumentError Base.delete(empty, 129)
end

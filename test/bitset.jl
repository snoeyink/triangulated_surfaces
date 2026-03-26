@testset "BitSet128" begin
    @test isbitstype(BitSet128)

    empty = BitSet128()
    @test empty.bits == zero(UInt128)
    @test iszero(empty)
    @test length(empty) == 0

    low = push(empty, 1)
    mid = push(empty, 64)
    high = push(empty, 128)
    pair = push(low, 128)
    bitxor = getproperty(Base, Symbol("⊻"))

    @test 1 in low
    @test !(2 in low)
    @test 64 in mid
    @test 128 in high
    @test length(BitSet128.add(low, 1)) == 1

    @test (low & high).bits == zero(UInt128)
    @test (low | high).bits == pair.bits
    @test bitxor(pair, low).bits == high.bits

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
    @test bitxor(a, b).bits == UInt128(0b0110)

    @test_throws ArgumentError (0 in empty)
    @test_throws ArgumentError (129 in empty)
    @test_throws ArgumentError push(empty, 0)
    @test_throws ArgumentError push(empty, 129)
    @test_throws ArgumentError delete(empty, 0)
    @test_throws ArgumentError delete(empty, 129)
end

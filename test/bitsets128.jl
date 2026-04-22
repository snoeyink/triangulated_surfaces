using Test
include(joinpath(@__DIR__, "..", "src", "BitSets128.jl"))

@testset "BitSets128 Module" begin

    @testset "BitSet128: Basics & Constructors" begin
        empty_set = BitSet128()
        @test isempty(empty_set)
        @test minimum(empty_set) == 0
        @test collect(empty_set) == Int[]

        s1 = singleton(BitSet128, 1)
        s64 = singleton(BitSet128, 64)
        s65 = singleton(BitSet128, 65)
        s128 = singleton(BitSet128, 128)

        @test 1 in s1
        @test 64 in s64
        @test 65 in s65
        @test 128 in s128
        @test !(2 in s1)
        
        @test minimum(s65) == 65
        @test !isempty(s1)
    end

    @testset "BitSet128: Bitwise Operations & Queries" begin
        s1 = singleton(BitSet128, 1)
        s2 = singleton(BitSet128, 2)
        s65 = singleton(BitSet128, 65)
        
        # Binary Ops
        u = s1 | s2 | s65
        @test 1 in u && 2 in u && 65 in u
        @test minimum(u) == 1
        
        @test (u & s2) == s2
        @test isempty(s1 & s2)
        @test isdisjoint(s1, s2)
        @test !isdisjoint(u, s1)
        
        @test setdiff(u, s1) == (s2 | s65)
        @test xor(u, s2) == (s1 | s65)
        
        # Unary Op
        inv_u = ~u
        @test !(1 in inv_u)
        @test 3 in inv_u
        @test 128 in inv_u
    end

    @testset "BitSet128: Valid Masking" begin
        m0 = valid_mask(BitSet128, 0)
        @test isempty(m0)
        
        m10 = valid_mask(BitSet128, 10)
        @test 10 in m10
        @test !(11 in m10)
        @test collect(m10) == collect(1:10)
        
        m64 = valid_mask(BitSet128, 64)
        @test 64 in m64
        @test !(65 in m64)
        
        m128 = valid_mask(BitSet128, 128)
        @test 128 in m128
        @test minimum(~m128) == 0 # ~m128 should be empty
        
        m_over = valid_mask(BitSet128, 200)
        @test m_over == m128
    end

    @testset "BitSet128: Iteration" begin
        s = singleton(BitSet128, 5) | singleton(BitSet128, 70) | singleton(BitSet128, 120)
        elements = collect(s)
        @test elements == [5, 70, 120]
    end

    @testset "BitSetOriented128: Basics & Iteration" begin
        empty_ori = BitSetOriented128()
        @test isempty(empty_ori)
        @test minimum(empty_ori) == 0
        
        # Create a mock oriented set manually via underlying BitSet128s
        fwd_set = singleton(BitSet128, 10) | singleton(BitSet128, 128)
        rev_set = singleton(BitSet128, 1)  | singleton(BitSet128, 127) # Represents 129 and 255
        ori = BitSetOriented128(fwd_set, rev_set)
        
        @test !isempty(ori)
        @test 10 in ori
        @test 128 in ori
        @test 129 in ori  # rev bit 1
        @test 255 in ori  # rev bit 127
        @test !(127 in ori)

        @test minimum(ori) == 10
        @test collect(ori) == [10, 128, 129, 255]
    end

    @testset "BitSetOriented128: Valid Masking" begin
        m10 = valid_mask(BitSetOriented128, 10)
        @test 10 in m10
        @test !(11 in m10)
        @test 138 in m10  # 128 + 10 in rev
        @test !(139 in m10)
        @test collect(m10) == [collect(1:10); collect(129:138)]
    end

    @testset "Rev Wrapper: Semantics" begin
        fwd1 = singleton(BitSet128, 5)
        rev1 = singleton(BitSet128, 10)
        o1 = BitSetOriented128(fwd1, rev1) # 5 and 138

        fwd2 = singleton(BitSet128, 10)
        rev2 = singleton(BitSet128, 5)
        o2 = BitSetOriented128(fwd2, rev2) # 10 and 133

        # Rev simply intercepts and pre-swaps o2's fwd/rev
        # So o2 "reversed" looks exactly like o1
        @test o1 == Rev(o2)
        @test Rev(o1) == o2
        @test Rev(o1) == Rev(o1)
        
        # Test Binary Ops with Rev
        u = o1 | Rev(o2)
        # Since Rev(o2) == o1, o1 | Rev(o2) should just be o1
        @test u == o1
        
        x = o1 & Rev(o1) 
        # o1 has fwd=5, rev=10. Rev(o1) acts as fwd=10, rev=5
        # Intersection should be empty
        @test isempty(x)
        @test isdisjoint(o1, Rev(o1))
        
        # Test Unary Op with Rev
        # ~Rev(o2) should be equivalent to ~(Rev(o2)) -> ~o1
        @test ~Rev(o2) == ~o1
    end

end
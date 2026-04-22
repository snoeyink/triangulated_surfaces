# test/bitsetoriented.jl
#
# Run standalone:  julia --project -e 'include("test/bitsetoriented.jl")'
# Or via Pkg.test() if included in test/runtests.jl.

using Test
include("../src/BitSetsOriented128.jl")
using .BitSetsOriented128

# ── helper: reconstruct a BitSetOriented128 from collected positions ──────────
function from_positions(iter)
    foldl(iter; init = BitSetOriented128()) do acc, pos
        if pos ≤ 128
            BitSetOriented128(acc.lo | one(UInt128) << (pos - 1), acc.hi)
        else
            BitSetOriented128(acc.lo, acc.hi | one(UInt128) << (pos - 129))
        end
    end
end

@testset "BitSetsOriented128" begin

    # ── edge_pos ──────────────────────────────────────────────────────────────
    @testset "edge_pos" begin
        # Manual evaluation of i + (j-1)(j-2)/2
        @test edge_pos(1, 2) === UInt8(1)
        @test edge_pos(1, 3) === UInt8(2)
        @test edge_pos(2, 3) === UInt8(3)
        @test edge_pos(1, 4) === UInt8(4)
        @test edge_pos(2, 4) === UInt8(5)
        @test edge_pos(3, 4) === UInt8(6)
        @test edge_pos(1, 5) === UInt8(7)
        # Maximum valid position for n=16 must fit in UInt8 and in UInt128
        @test edge_pos(15, 16) === UInt8(120)
        # Return type contract
        @test edge_pos(1, 2) isa UInt8
    end

    # ── singleton ─────────────────────────────────────────────────────────────
    @testset "singleton" begin
        # Forward edge 1→2: lo bit at position 1 (index 0), hi empty
        s12 = singleton(1, 2)
        @test s12.lo === one(UInt128) << 0
        @test s12.hi === zero(UInt128)

        # Reverse edge 2→1: same bit index but in hi
        s21 = singleton(2, 1)
        @test s21.lo === zero(UInt128)
        @test s21.hi === one(UInt128) << 0

        # Both directions share the same bit position
        @test singleton(1, 2).lo === singleton(2, 1).hi
        @test singleton(3, 7).lo === singleton(7, 3).hi

        # Forward edge 2→3: edge_pos(2,3) = 3 → bit index 2
        s23 = singleton(2, 3)
        @test s23.lo === one(UInt128) << 2
        @test s23.hi === zero(UInt128)

        # Maximum valid edge: 15→16, edge_pos = 120 → bit index 119
        s_max = singleton(15, 16)
        @test s_max.lo === one(UInt128) << 119
        @test s_max.hi === zero(UInt128)
        @test singleton(16, 15).hi === one(UInt128) << 119

        # Bounds errors
        @test_throws ArgumentError singleton(0, 1)     # vertex < 1
        @test_throws ArgumentError singleton(1, 17)    # vertex > 16
        @test_throws ArgumentError singleton(16, 17)   # vertex > 16
        @test_throws ArgumentError singleton(1, 1)     # self-loop
    end

    # ── constructors and zero ─────────────────────────────────────────────────
    @testset "constructors" begin
        @test iszero(BitSetOriented128())
        @test iszero(zero(BitSetOriented128))
        @test BitSetOriented128() == zero(BitSetOriented128)
        @test BitSetOriented128().lo === zero(UInt128)
        @test BitSetOriented128().hi === zero(UInt128)
    end

    # ── triangle ──────────────────────────────────────────────────────────────
    @testset "triangle" begin
        # triangle(1,2,3): 1→2 (lo pos 1), 2→3 (lo pos 3), 3→1 (hi pos 2)
        #   lo: bit indices 0 and 2 → 0b101 = 5
        #   hi: bit index  1       → 0b010 = 2
        t = triangle(1, 2, 3)
        @test t.lo === UInt128(0b101)
        @test t.hi === UInt128(0b010)
        @test length(t) == 3

        # Reversed triangle swaps lo ↔ hi
        @test reversed(t).lo === t.hi
        @test reversed(t).hi === t.lo

        # A triangle and its reverse are complementary within valid_mask(3)
        @test isdisjoint(t, reversed(t))
        @test (t | reversed(t)) == valid_mask(3)

        # Cyclic symmetry: all three rotations are equal
        @test triangle(1, 2, 3) == triangle(2, 3, 1)
        @test triangle(1, 2, 3) == triangle(3, 1, 2)

        # Opposite winding is not equal
        @test triangle(1, 2, 3) != triangle(1, 3, 2)

        # Opposite winding equals reversed
        @test triangle(1, 3, 2) == reversed(triangle(1, 2, 3))
    end

    # ── reversed and rev ──────────────────────────────────────────────────────
    @testset "reversed and rev" begin
        s = singleton(1, 3) | singleton(3, 2)   # one forward, one reverse

        # reversed materialises a new set with lo ↔ hi
        r = reversed(s)
        @test r.lo === s.hi
        @test r.hi === s.lo

        # reversed is an involution
        @test reversed(reversed(s)) == s

        # reversed of empty is empty
        @test reversed(BitSetOriented128()) == BitSetOriented128()

        # rev wraps lazily (same object inside)
        @test rev(s).s === s

        # rev in binary ops is identical to materialised reversed
        a = triangle(1, 2, 3)
        b = triangle(2, 3, 4)
        @test (a ⊻ rev(b)) == (a ⊻ reversed(b))
        @test (a &  rev(b)) == (a &  reversed(b))
        @test (a |  rev(b)) == (a |  reversed(b))

        # ~rev(b): lazy complement-of-reverse
        #   ~Rev(b) == Rev(~b), so (a & ~rev(b)) == (a & reversed(~b))
        @test (~rev(b)).s == ~b
        @test (a & ~rev(b)) == (a & reversed(~b))
    end

    # ── bitwise operations ────────────────────────────────────────────────────
    @testset "bitwise" begin
        s12 = singleton(1, 2)
        s21 = singleton(2, 1)
        s13 = singleton(1, 3)

        # OR
        u = s12 | s21
        @test u.lo === s12.lo
        @test u.hi === s21.hi
        @test length(u) == 2

        # AND
        @test (s12 & (s12 | s13)) == s12
        @test (s12 & s21)          == BitSetOriented128()   # lo vs hi, no overlap

        # XOR
        @test ((s12 | s21) ⊻ s12) == s21
        @test (s12 ⊻ s12)          == BitSetOriented128()

        # NOT: bits 121–128 are collateral damage; masking cleans them
        @test (~BitSetOriented128() & valid_mask(3)) == valid_mask(3)

        # Complement within valid_mask is a true set complement
        t   = triangle(1, 2, 3)
        cmp = ~t & valid_mask(3)
        @test isdisjoint(t, cmp)
        @test (t | cmp) == valid_mask(3)
    end

    # ── predicates ────────────────────────────────────────────────────────────
    @testset "predicates" begin
        empty_s = BitSetOriented128()
        s12     = singleton(1, 2)
        s21     = singleton(2, 1)

        # iszero / isempty
        @test  iszero(empty_s)
        @test  isempty(empty_s)
        @test !iszero(s12)
        @test !isempty(s12)

        # length
        @test length(empty_s)       == 0
        @test length(s12)           == 1
        @test length(s12 | s21)     == 2
        @test length(triangle(1,2,3)) == 3

        # issubset
        @test  issubset(empty_s, s12)
        @test  issubset(empty_s, empty_s)
        @test  issubset(s12, s12)
        @test  issubset(s12, s12 | s21)
        @test !issubset(s12 | s21, s12)

        # isdisjoint
        @test  isdisjoint(empty_s, s12)
        @test  isdisjoint(s12, s21)           # forward vs reverse, different fields
        @test  isdisjoint(s12, singleton(1, 3))
        @test !isdisjoint(s12, s12 | s21)
        @test  isdisjoint(triangle(1,2,3), reversed(triangle(1,2,3)))
    end

    # ── min_pos ───────────────────────────────────────────────────────────────
    @testset "min_pos" begin
        # Empty: sentinel 0
        @test min_pos(BitSetOriented128()) === UInt16(0)

        # lo bits: result == edge_pos
        @test min_pos(singleton(1, 2)) === UInt16(1)    # edge_pos(1,2) = 1
        @test min_pos(singleton(2, 3)) === UInt16(3)    # edge_pos(2,3) = 3

        # hi bits only: result == 128 + edge_pos
        @test min_pos(singleton(2, 1)) === UInt16(129)  # 128 + edge_pos(1,2) = 129
        @test min_pos(singleton(3, 2)) === UInt16(131)  # 128 + edge_pos(2,3) = 131

        # Maximum edge: 15→16 in lo, 16→15 in hi
        @test min_pos(singleton(15, 16)) === UInt16(120)
        @test min_pos(singleton(16, 15)) === UInt16(248)  # 128 + 120

        # lo is always checked before hi, regardless of bit value
        # lo at pos 3, hi at pos 1 → lo wins
        mixed = singleton(2, 3) | singleton(2, 1)
        @test min_pos(mixed) === UInt16(3)

        # Multiple lo bits: lowest wins
        @test min_pos(singleton(2,3) | singleton(1,2)) === UInt16(1)

        # Multiple hi bits: lowest wins (both hi; min edge_pos is 1 → result 129)
        @test min_pos(singleton(3,1) | singleton(2,1)) === UInt16(129)

        # Return type
        @test min_pos(singleton(1, 2)) isa UInt16
    end

    # ── set_both ──────────────────────────────────────────────────────────────
    @testset "set_both" begin
        empty_s = BitSetOriented128()

        # set_both(s, pos): sets the same bit in both lo and hi
        s1 = set_both(empty_s, 1)
        @test s1.lo === UInt128(1)
        @test s1.hi === UInt128(1)

        # set_both(s, a, b): convenience form; edge_pos(1,2) = 1
        s2 = set_both(empty_s, 1, 2)
        @test s2 == s1

        # Unrelated bits in lo are preserved; hi bit is newly set
        base = singleton(1, 3)                  # lo bit at pos 2; hi empty
        s3   = set_both(base, 1, 2)
        @test s3.lo === base.lo | UInt128(1)    # bit at pos 1 added
        @test s3.hi === UInt128(1)

        # Idempotent
        @test set_both(s1, 1)    == s1
        @test set_both(s2, 1, 2) == s2

        # Bounds errors
        @test_throws ArgumentError set_both(empty_s, 0, 1)
        @test_throws ArgumentError set_both(empty_s, 1, 17)
    end

    # ── in_directed, push, delete ─────────────────────────────────────────────
    @testset "membership, push, delete" begin
        t = triangle(1, 2, 3)   # 1→2, 2→3, 3→1

        # in_directed: present edges
        @test  in_directed(1, 2, t)
        @test  in_directed(2, 3, t)
        @test  in_directed(3, 1, t)

        # in_directed: absent (reverse) edges
        @test !in_directed(2, 1, t)
        @test !in_directed(3, 2, t)
        @test !in_directed(1, 3, t)

        # push adds a directed edge
        t2 = push(t, 2, 1)
        @test  in_directed(2, 1, t2)
        @test  in_directed(1, 2, t2)   # original preserved
        @test  length(t2) == 4

        # push is idempotent
        @test push(t, 1, 2) == t

        # delete removes a directed edge
        t3 = delete(t, 1, 2)
        @test !in_directed(1, 2, t3)
        @test  in_directed(2, 3, t3)   # other edges unaffected
        @test  in_directed(3, 1, t3)
        @test  length(t3) == 2

        # delete of absent edge is identity
        @test delete(t, 2, 1) == t
    end

    # ── iteration ─────────────────────────────────────────────────────────────
    @testset "iteration" begin
        # eltype contract
        @test eltype(BitSetOriented128) == UInt16

        # Empty set
        @test collect(BitSetOriented128()) == UInt16[]

        # Single forward edge
        @test collect(singleton(1, 2)) == [UInt16(1)]

        # Single reverse edge: encoded as 128 + edge_pos
        @test collect(singleton(2, 1)) == [UInt16(129)]

        # triangle(1,2,3): lo at {pos 1, pos 3}, hi at {pos 2 → 130}
        t = triangle(1, 2, 3)
        @test collect(t) == [UInt16(1), UInt16(3), UInt16(130)]

        # Length and collect agree
        big = valid_mask(5)   # C(5,2)=10 undirected → 20 directed
        @test length(big) == 20
        @test length(collect(big)) == 20

        # Positions are strictly increasing (lo ascending, then hi ascending)
        @test issorted(collect(valid_mask(4)))

        # Round-trip: iterating and reconstructing gives back the original
        @test from_positions(t)              == t
        @test from_positions(valid_mask(4))  == valid_mask(4)
        @test from_positions(BitSetOriented128()) == BitSetOriented128()
    end

    # ── valid_mask ────────────────────────────────────────────────────────────
    @testset "valid_mask" begin
        @test valid_mask(0) == BitSetOriented128()

        # n=1: no edges
        @test valid_mask(1) == BitSetOriented128()

        # n=2: 1 undirected edge → 2 directed
        vm2 = valid_mask(2)
        @test vm2.lo === UInt128(1)
        @test vm2.hi === UInt128(1)
        @test length(vm2) == 2

        # n=3: 3 undirected → 6 directed
        vm3 = valid_mask(3)
        @test vm3.lo === UInt128(0b111)
        @test vm3.hi === UInt128(0b111)
        @test length(vm3) == 6

        # n=4: C(4,2)=6 undirected → 12 directed
        @test length(valid_mask(4)) == 12

        # lo == hi for all n (symmetric position scheme)
        for n in 0:16
            vm = valid_mask(n)
            @test vm.lo == vm.hi
        end

        # n=16: C(16,2)=120 undirected → 240 directed
        vm16 = valid_mask(16)
        @test count_ones(vm16.lo) == 120
        @test count_ones(vm16.hi) == 120
        @test length(vm16) == 240

        # Masked complement is a true complement
        t   = triangle(1, 2, 3)
        cmp = ~t & valid_mask(3)
        @test isdisjoint(t, cmp)
        @test (t | cmp) == vm3

        # Bounds errors
        @test_throws ArgumentError valid_mask(-1)
        @test_throws ArgumentError valid_mask(17)
    end

    # ── xor_rev / ⊻ rev semantics ─────────────────────────────────────────────
    @testset "xor_rev" begin
        t = triangle(1, 2, 3)

        # t and reversed(t) are disjoint, so XOR == OR here
        @test (t ⊻ rev(t)) == (t | reversed(t))

        # A set equal to its own reverse XORed with itself-reversed gives zero
        both = singleton(1, 2) | singleton(2, 1)   # lo=1, hi=1; reversed == itself
        @test iszero(both ⊻ rev(both))

        # Lazy rev matches eager reversed in all three operations
        a = triangle(1, 2, 4)
        b = triangle(2, 3, 4)
        @test (a ⊻ rev(b)) == (a ⊻ reversed(b))
        @test (a &  rev(b)) == (a &  reversed(b))
        @test (a |  rev(b)) == (a |  reversed(b))

        # a & ~rev(b): keep edges of a whose reverse is absent from b.
        # t and reversed(t) are disjoint, so ~reversed(t) ⊇ t, giving back t.
        @test (t & ~rev(t)) == t

        # Removing edges of a that also appear reversed in b
        a2 = singleton(1,2) | singleton(2,3)   # forward edges
        b2 = singleton(2,1)                    # reverse of 1→2
        # a2 & ~rev(b2) should drop 1→2 (its reverse 2→1 is in b2)
        @test (a2 & ~rev(b2)) == singleton(2, 3)
    end

    # ── integration: tournament on 4 vertices ─────────────────────────────────
    @testset "tournament K4" begin
        # Transitive tournament: i→j for all 1 ≤ i < j ≤ 4
        tour = foldl(|, [singleton(i, j) for i in 1:4 for j in (i+1):4])

        @test length(tour) == 6            # C(4,2)

        # Tournament ∪ reversed tournament = all directed edges on 4 vertices
        @test (tour | reversed(tour)) == valid_mask(4)

        # No reciprocal edges
        @test isdisjoint(tour, reversed(tour))

        # Removing and re-adding one edge is an identity
        e = singleton(2, 4)
        @test push(delete(tour, 2, 4), 2, 4) == tour

        # Iterating gives exactly 6 positions, all in lo (all forward edges)
        positions = collect(tour)
        @test length(positions) == 6
        @test all(p -> p ≤ 128, positions)   # no hi bits in a transitive tournament
        @test issorted(positions)
    end

end
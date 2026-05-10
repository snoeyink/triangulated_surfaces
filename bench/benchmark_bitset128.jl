# benchmark_bitset128.jl
#
# Compare: 2 × UInt128  vs  4 × UInt64  for oriented-edge bitsets.
# Test:    pairwise union (|) over 2^12 random sets,
#          each with 3–6 bits set at distinct positions in 1..120.

using BenchmarkTools
using Random

# ─────────────────────────────────────────────────────────────────────────────
# Implementation A: 2 × UInt128
# ─────────────────────────────────────────────────────────────────────────────

struct OrientedSet2x128
    lo::UInt128
    hi::UInt128
end

@inline Base.:|(a::OrientedSet2x128, b::OrientedSet2x128) =
    OrientedSet2x128(a.lo | b.lo, a.hi | b.hi)

# ─────────────────────────────────────────────────────────────────────────────
# Implementation B: 4 × UInt64
# (lo split into lo_lo / lo_hi, same for hi)
# ─────────────────────────────────────────────────────────────────────────────

struct OrientedSet4x64
    lo_lo::UInt64        # bits  0–63  of the forward-edge word
    lo_hi::UInt64        # bits 64–119 of the forward-edge word
    hi_lo::UInt64        # bits  0–63  of the reverse-edge word
    hi_hi::UInt64        # bits 64–119 of the reverse-edge word
end

@inline Base.:|(a::OrientedSet4x64, b::OrientedSet4x64) =
    OrientedSet4x64(
        a.lo_lo | b.lo_lo,
        a.lo_hi | b.lo_hi,
        a.hi_lo | b.hi_lo,
        a.hi_hi | b.hi_hi,
    )

# ─────────────────────────────────────────────────────────────────────────────
# Data generation
# ─────────────────────────────────────────────────────────────────────────────

"""
    sparse_bits(rng, n) -> UInt128

UInt128 with exactly `n` bits set at distinct positions drawn uniformly
from 0..119 (the 120 undirected-edge slots for ≤ 16 vertices).
"""
function sparse_bits(rng::AbstractRNG, n::Int)::UInt128
    bits  = zero(UInt128)
    count = 0
    while count < n
        pos  = rand(rng, 0:119)
        mask = one(UInt128) << pos
        if iszero(bits & mask)
            bits  |= mask
            count += 1
        end
    end
    bits
end

const N = 1 << 12   # 4 096

function make_sets_2x128(rng::AbstractRNG)::Vector{OrientedSet2x128}
    sets = Vector{OrientedSet2x128}(undef, N)
    for i in 1:N
        sets[i] = OrientedSet2x128(
            sparse_bits(rng, rand(rng, 3:6)),
            sparse_bits(rng, rand(rng, 3:6)),
        )
    end
    sets
end

"""Convert to 4×UInt64 layout from the same underlying data (fair comparison)."""
function to_4x64(sets::Vector{OrientedSet2x128})::Vector{OrientedSet4x64}
    map(sets) do s
        OrientedSet4x64(
            s.lo % UInt64,           # lower 64 bits of lo
            (s.lo >> 64) % UInt64,   # upper 64 bits of lo
            s.hi % UInt64,           # lower 64 bits of hi
            (s.hi >> 64) % UInt64,   # upper 64 bits of hi
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Benchmark kernel  (generic over both types — one function, two specialisations)
# ─────────────────────────────────────────────────────────────────────────────

"""
    union_pairs!(out, sets) -> out

Fill `out[k]` with `sets[2k-1] | sets[2k]` for k = 1 .. N÷2.
"""
function union_pairs!(out::Vector{T}, sets::Vector{T}) where {T}
    @inbounds for i in eachindex(sets), j in eachindex(sets)
        out[i] = sets[i] | sets[j]
    end
    out
end

# ─────────────────────────────────────────────────────────────────────────────
# Setup  (runs once at include-time)
# ─────────────────────────────────────────────────────────────────────────────

const rng        = MersenneTwister(42)
const sets_2x128 = make_sets_2x128(rng)
const sets_4x64  = to_4x64(sets_2x128)          # identical data, different layout

const out_2x128  = Vector{OrientedSet2x128}(undef, N )
const out_4x64   = Vector{OrientedSet4x64}(undef, N)

# ── Zero-allocation guard ─────────────────────────────────────────────────────
let a = @allocated(union_pairs!(out_2x128, sets_2x128)),
    b = @allocated(union_pairs!(out_4x64,  sets_4x64))
    @assert a == 0 "2×UInt128 kernel allocated $a bytes — check type stability"
    @assert b == 0 "4×UInt64  kernel allocated $b bytes — check type stability"
end

# ── Correctness check ─────────────────────────────────────────────────────────
@assert all(eachindex(out_2x128)) do k
    s = out_2x128[k];  t = out_4x64[k]
    s.lo % UInt64        == t.lo_lo  &&
    (s.lo >> 64) % UInt64 == t.lo_hi  &&
    s.hi % UInt64        == t.hi_lo  &&
    (s.hi >> 64) % UInt64 == t.hi_hi
end "Implementations disagree — bug in conversion or kernel"

# ─────────────────────────────────────────────────────────────────────────────
# Benchmarks
# ─────────────────────────────────────────────────────────────────────────────

println("\n=== 2 × UInt128  ($(N^2) unions) ===")
display(@benchmark union_pairs!($out_2x128, $sets_2x128))

println("\n=== 4 × UInt64   ($(N^2) unions) ===")
display(@benchmark union_pairs!($out_4x64,  $sets_4x64))

# For full statistics (min / median / mean / allocs / GC), swap @btime for:
#   display(@benchmark union_pairs!($out_2x128, $sets_2x128))
#   display(@benchmark union_pairs!($out_4x64,  $sets_4x64))
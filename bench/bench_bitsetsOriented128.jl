# bench_bitsetsOriented128.jl
#
# Compare BitSetOriented128 implementations:
#   A) from src/BitSets128.jl
#   B) from src/BitSetsOriented128.jl
#
# Benchmark kernel: N^2 unions over N random sparse sets.

using BenchmarkTools
using Random

include("../src/BitSets128.jl")
include("../src/BitSetsOriented128.jl")

const BS128 = Main
const BSO128 = Main.BitSetsOriented128

# Data generation helpers

"""
    sparse_bits(rng, n) -> UInt128

UInt128 with exactly `n` bits set at distinct positions in 1..120.
"""
function sparse_bits(rng::AbstractRNG, n::Int)::UInt128
    bits = zero(UInt128)
    count = 0
    while count < n
        pos = rand(rng, 1:120)
        mask = one(UInt128) << (pos - 1)
        if iszero(bits & mask)
            bits |= mask
            count += 1
        end
    end
    return bits
end

"""
    make_masks(rng, nsets) -> Vector{Tuple{UInt128,UInt128}}

Creates `(lo, hi)` masks where each side has 3-6 bits set in positions 1..120.
"""
function make_masks(rng::AbstractRNG, nsets::Int)::Vector{Tuple{UInt128,UInt128}}
    masks = Vector{Tuple{UInt128,UInt128}}(undef, nsets)
    for i in 1:nsets
        lo = sparse_bits(rng, rand(rng, 3:6))
        hi = sparse_bits(rng, rand(rng, 3:6))
        masks[i] = (lo, hi)
    end
    return masks
end

# Build implementation A sets from shared masks
function to_sets_bs128(masks::Vector{Tuple{UInt128,UInt128}})::Vector{BS128.BitSetOriented128}
    sets = Vector{BS128.BitSetOriented128}(undef, length(masks))
    for i in eachindex(masks)
        lo, hi = masks[i]
        s = BS128.BitSetOriented128()
        for p in 1:120
            bit = one(UInt128) << (p - 1)
            if !iszero(lo & bit)
                s |= BS128.singleton(BS128.BitSetOriented128, p)
            end
            if !iszero(hi & bit)
                s |= BS128.singleton(BS128.BitSetOriented128, p + 128)
            end
        end
        sets[i] = s
    end
    return sets
end

# Build implementation B sets from shared masks
function to_sets_bso128(masks::Vector{Tuple{UInt128,UInt128}})::Vector{BSO128.BitSetOriented128}
    sets = Vector{BSO128.BitSetOriented128}(undef, length(masks))
    for i in eachindex(masks)
        lo, hi = masks[i]
        sets[i] = BSO128.BitSetOriented128(lo, hi)
    end
    return sets
end

# Convert implementation A to (lo, hi) masks for cross-checking
@inline function split_bs128(s::BS128.BitSetOriented128)::Tuple{UInt128,UInt128}
    lo = UInt128(s.fwd.words[1]) | (UInt128(s.fwd.words[2]) << 64)
    hi = UInt128(s.rev.words[1]) | (UInt128(s.rev.words[2]) << 64)
    return (lo, hi)
end

# Benchmark kernel (same algorithm for both implementations)
"""
    union_allpairs!(out, sets) -> out

Fill `out[i]` repeatedly with `sets[i] | sets[j]` for all `i, j`.
"""
function union_allpairs!(out::Vector{T}, sets::Vector{T}) where {T}
    @inbounds for i in eachindex(sets), j in eachindex(sets)
        out[i] = sets[i] | sets[j]
    end
    return out
end

# Setup
const N = 1 << 12  # 4_096
const rng = MersenneTwister(42)
const masks = make_masks(rng, N)

const sets_bs128 = to_sets_bs128(masks)
const sets_bso128 = to_sets_bso128(masks)

const out_bs128 = Vector{BS128.BitSetOriented128}(undef, N)
const out_bso128 = Vector{BSO128.BitSetOriented128}(undef, N)

# Zero-allocation guard
let a = @allocated(union_allpairs!(out_bs128, sets_bs128)),
    b = @allocated(union_allpairs!(out_bso128, sets_bso128))
    @assert a == 0 "BitSets128.BitSetOriented128 kernel allocated $a bytes"
    @assert b == 0 "BitSetsOriented128.BitSetOriented128 kernel allocated $b bytes"
end

# Correctness check
@assert all(eachindex(out_bs128)) do k
    split_bs128(out_bs128[k]) == (out_bso128[k].lo, out_bso128[k].hi)
end "Implementations disagree on union_allpairs!"

# Benchmarks
println("\n=== BitSets128.BitSetOriented128  ($(N^2) unions) ===")
display(@benchmark union_allpairs!($out_bs128, $sets_bs128))

println("\n=== BitSetsOriented128.BitSetOriented128  ($(N^2) unions) ===")
display(@benchmark union_allpairs!($out_bso128, $sets_bso128))

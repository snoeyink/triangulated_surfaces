using BenchmarkTools
using Random

include(joinpath(@__DIR__, "..", "src", "BitSets2n.jl"))
include(joinpath(@__DIR__, "..", "src", "BitSets128.jl"))

const BS2 = BitSets2n
const BS64 = BitSets128

const N = 1024
const SEED = 42

@inline function as_bs64(bits::UInt128)
    lo = UInt64(bits & UInt128(typemax(UInt64)))
    hi = UInt64(bits >> 64)
    BS64.BitSet128((lo, hi))
end

function make_inputs(n::Int, seed::Int)
    rng = MersenneTwister(seed)
    raw_a = rand(rng, UInt128, n)
    raw_b = rand(rng, UInt128, n)
    elems = rand(rng, 1:128, n)

    a2 = [BS2.BitSet128(x) for x in raw_a]
    b2 = [BS2.BitSet128(x) for x in raw_b]
    a64 = [as_bs64(x) for x in raw_a]
    b64 = [as_bs64(x) for x in raw_b]

    return a2, b2, a64, b64, elems
end

@inline setdiff_like(a, b) = a & ~b

function k_or(a, b)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = a[i] | b[i]
    end
    x
end

function k_and(a, b)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = a[i] & b[i]
    end
    x
end

function k_xor_2n(a, b)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = a[i] ⊻ b[i]
    end
    x
end

function k_xor_64(a, b)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = xor(a[i], b[i])
    end
    x
end

function k_not(a)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = ~a[i]
    end
    x
end

function k_setdiff(a, b)
    x = a[1]
    @inbounds for i in eachindex(a)
        x = setdiff_like(a[i], b[i])
    end
    x
end

function k_eq(a, b)
    s = 0
    @inbounds for i in eachindex(a)
        s += Int(a[i] == b[i])
    end
    s
end

function k_isempty(a)
    s = 0
    @inbounds for i in eachindex(a)
        s += Int(isempty(a[i]))
    end
    s
end

function k_isdisjoint(a, b)
    s = 0
    @inbounds for i in eachindex(a)
        s += Int(isdisjoint(a[i], b[i]))
    end
    s
end

function k_in(a, elems)
    s = 0
    @inbounds for i in eachindex(a)
        s += Int(in(elems[i], a[i]))
    end
    s
end

function k_iter(a)
    s = 0
    @inbounds for i in eachindex(a)
        for x in a[i]
            s += Int(x)
        end
    end
    s
end

function run_benchmarks()
    a2, b2, a64, b64, elems = make_inputs(N, SEED)

    println("BitSet128 benchmark (N=$(N), seed=$(SEED))")
    println("Comparing BitSets2n.BitSet128 (UInt128) vs BitSets128.BitSet128 (2xUInt64)")

    println("\n=== BitSets2n (UInt128-backed) ===")
    display(@benchmark k_or($a2, $b2))
    display(@benchmark k_and($a2, $b2))
    display(@benchmark k_xor_2n($a2, $b2))
    display(@benchmark k_not($a2))
    display(@benchmark k_setdiff($a2, $b2))
    display(@benchmark k_eq($a2, $b2))
    display(@benchmark k_isempty($a2))
    display(@benchmark k_isdisjoint($a2, $b2))
    display(@benchmark k_in($a2, $elems))
    display(@benchmark k_iter($a2))

    println("\n=== BitSets128 (UInt64 tuple-backed) ===")
    display(@benchmark k_or($a64, $b64))
    display(@benchmark k_and($a64, $b64))
    display(@benchmark k_xor_64($a64, $b64))
    display(@benchmark k_not($a64))
    display(@benchmark k_setdiff($a64, $b64))
    display(@benchmark k_eq($a64, $b64))
    display(@benchmark k_isempty($a64))
    display(@benchmark k_isdisjoint($a64, $b64))
    display(@benchmark k_in($a64, $elems))
    display(@benchmark k_iter($a64))

    nothing
end

run_benchmarks()

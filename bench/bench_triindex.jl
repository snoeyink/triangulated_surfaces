# bench/bench_triindex.jl
# Compares 3-D array indexing strategies for tri_table access.
# Four cases: n×n×n Array, 16×16×16 Array,
#             flat n³ Vector (const stride), flat 16³ Vector (bit-shifts).

using BenchmarkTools
using Random

const n = 13    # vertex count under test
const P = 16    # power-of-2 upper bound (MAX_VERTICES)

# ─── Kernels ──────────────────────────────────────────────────────────────────

# Cases 1 & 2: Array — strides are runtime values; LLVM sees generic multiply.
# Separate functions so @benchmark can specialize, but generated code is identical.
function k_array_n!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += Int16(1)
    end
end

function k_array_16!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += Int16(1)
    end
end

# Case 3: Flat n³ vector; stride n=13 is a compile-time constant.
# LLVM can strength-reduce: n*x = (x<<3) + (x<<2) + x, etc.
function k_flat_n!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i] - 1) + n * (r[i+1] - 1) + n * n * (r[i+2] - 1)
        tri[idx + 1] += Int16(1)
    end
end

# Case 4: Flat 16³ vector; index built with pure bit-ops (no multiply).
# Valid because values in [1..13] ⊂ [1..16], so (x-1) fits in 4 bits.
# Layout: bit 0-3 = first dim, bits 4-7 = second, bits 8-11 = third.
function k_flat_16!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i] - 1) | ((r[i+1] - 1) << 4) | ((r[i+2] - 1) << 8)
        tri[idx + 1] += Int16(1)
    end
end

# ─── Setup ────────────────────────────────────────────────────────────────────
Random.seed!(42)
const r_bench = rand(1:n, 1000)

arr_n    = zeros(Int16, n, n, n)     # 2197 elements, 4.3 KB
arr_16   = zeros(Int16, P, P, P)     # 4096 elements, 8.0 KB
flat_n   = zeros(Int16, n^3)         # 2197 elements, 4.3 KB
flat_16  = zeros(Int16, P^3)         # 4096 elements, 8.0 KB

# ─── Sanity check (run once before benchmarking) ──────────────────────────────
k_array_n!(arr_n, r_bench);   k_array_16!(arr_16, r_bench)
k_flat_n!(flat_n, r_bench);   k_flat_16!(flat_16, r_bench)

# Maxima must agree (same random data, same index mapping)
@assert maximum(arr_n)   == maximum(flat_n)  "flat_n indexing mismatch"
@assert maximum(arr_16)  == maximum(flat_16) "flat_16 indexing mismatch"
println("Sanity check passed.  max = $(maximum(arr_n))\n")

# ─── Benchmarks ───────────────────────────────────────────────────────────────
println("=== 1. Array{Int16,3} $(n)×$(n)×$(n)  (runtime stride=$(n)) ===")
display(@benchmark k_array_n!($arr_n, $r_bench))

println("\n=== 2. Array{Int16,3} $(P)×$(P)×$(P)  (runtime stride=$(P)) ===")
display(@benchmark k_array_16!($arr_16, $r_bench))

println("\n=== 3. flat Vector{Int16} length $(n^3), hardcoded stride=$n ===")
display(@benchmark k_flat_n!($flat_n, $r_bench))

println("\n=== 4. flat Vector{Int16} length $(P^3), bit-shift index ===")
display(@benchmark k_flat_16!($flat_16, $r_bench))

# ─── Optional: inspect generated LLVM ────────────────────────────────────────
# Uncomment to compare multiply vs shift in the inner loop:
# println("\n--- LLVM: k_flat_n! ---");  @code_llvm k_flat_n!(flat_n, r_bench)
# println("\n--- LLVM: k_flat_16! ---"); @code_llvm k_flat_16!(flat_16, r_bench)

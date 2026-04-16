# bench/bench_triindex.jl
#
# Six 3-D indexing strategies for tri_table.
# Write benchmark: 1000 × 998 increments, then find max.
# Read  benchmark: 1000 × 998 reads accumulated into a sum (closer to
#                  actual backtracking use; no mutation, no reset needed).
#
# Usage (VSCode):
#   1. Open integrated terminal: Ctrl+`
#   2. julia --project=.
#   3. include("bench/bench_triindex.jl")
#
# WARNING: Set RUN_MARRAY = false on a first run.
#          MArray{Tuple{13,13,13}} (2197 elements) and
#          MArray{Tuple{16,16,16}} (4096 elements) can take several minutes
#          to specialise; StaticArrays is designed for arrays < ~100 elements.

using BenchmarkTools
using Random
using StaticArrays

const RUN_MARRAY = true   # ← set false to skip cases 5 & 6

const n = 13
const P = 16

Random.seed!(42)
const r_bench = rand(1:n, 1000)   # values in 1..13; fits all six array shapes

# ─── Helper functions for MArray setup ────────────────────────────────────────
# Defined before kernels so @benchmark setup= can call them.
@inline fresh_mn()  = (m = MArray{Tuple{n,n,n}, UInt16}(undef); fill!(m, UInt16(0)); m)
@inline fresh_m16() = (m = MArray{Tuple{P,P,P}, UInt16}(undef); fill!(m, UInt16(0)); m)

# ─── Write kernels ─────────────────────────────────────────────────────────────

# 1: Array n×n×n  — stride loaded from array header; LLVM uses runtime multiply
function kw_array_n!(tri::Array{UInt16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += UInt16(1)
    end
end

# 2: Array 16×16×16 — same codegen path as case 1; different runtime stride
function kw_array_16!(tri::Array{UInt16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += UInt16(1)
    end
end

# 3: flat Vector, stride n=13 is a compile-time constant.
#    LLVM strength-reduces: 13x = (8+4+1)x → shifts+adds, no imul.
function kw_flat_n!(tri::Vector{UInt16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) + n*(r[i+1]-1) + n*n*(r[i+2]-1)
        tri[idx+1] += UInt16(1)
    end
end

# 4: flat Vector 16³, pure bit-shift index — zero multiplies.
#    Valid because (x-1) ∈ 0..15 for x ∈ 1..16, so each field fits in 4 bits.
#    Equivalent to column-major Array{16,16,16} indexing (verified in sanity check).
function kw_flat_16!(tri::Vector{UInt16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) | ((r[i+1]-1) << 4) | ((r[i+2]-1) << 8)
        tri[idx+1] += UInt16(1)
    end
end

# 5: MArray n×n×n — strides encoded in type parameters; LLVM sees constants.
#    Heap-allocated at this size, but no runtime stride lookup from header.
if RUN_MARRAY
    function kw_marray_n!(tri::MArray{Tuple{n,n,n}, UInt16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            tri[r[i], r[i+1], r[i+2]] += UInt16(1)
        end
    end

    # 6: MArray 16×16×16 — power-of-2 strides, compile-time
    function kw_marray_16!(tri::MArray{Tuple{P,P,P}, UInt16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            tri[r[i], r[i+1], r[i+2]] += UInt16(1)
        end
    end
end

# ─── Read kernels ──────────────────────────────────────────────────────────────
# Returns accumulated sum (prevents dead-code elimination; matches backtracking
# use pattern where tri_table is written once and read millions of times).

function kr_array_n(tri::Array{UInt16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        s += tri[r[i], r[i+1], r[i+2]]
    end
    return s
end

function kr_array_16(tri::Array{UInt16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        s += tri[r[i], r[i+1], r[i+2]]
    end
    return s
end

function kr_flat_n(tri::Vector{UInt16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) + n*(r[i+1]-1) + n*n*(r[i+2]-1)
        s += tri[idx+1]
    end
    return s
end

function kr_flat_16(tri::Vector{UInt16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) | ((r[i+1]-1) << 4) | ((r[i+2]-1) << 8)
        s += tri[idx+1]
    end
    return s
end

if RUN_MARRAY
    function kr_marray_n(tri::MArray{Tuple{n,n,n}, UInt16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            s += tri[r[i], r[i+1], r[i+2]]
        end
        return s
    end

    function kr_marray_16(tri::MArray{Tuple{P,P,P}, UInt16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            s += tri[r[i], r[i+1], r[i+2]]
        end
        return s
    end
end

# ─── Sanity checks ─────────────────────────────────────────────────────────────
println("Running sanity checks...")

# n-wide variants must agree on element values
let an = zeros(UInt16, n, n, n),
    fn = zeros(UInt16, n^3)
    kw_array_n!(an, r_bench);  kw_flat_n!(fn, r_bench)
    @assert vec(an) == vec(fn) "flat_n index mapping differs from Array n³"
    println("  Array n³ vs flat_n:  ✓  max=$(maximum(an))  sum=$(sum(an))")
end

# 16-wide variants must agree on element values
let a16 = zeros(UInt16, P, P, P),
    f16 = zeros(UInt16, P^3)
    kw_array_16!(a16, r_bench);  kw_flat_16!(f16, r_bench)
    @assert vec(a16) == vec(f16) "flat_16 index mapping differs from Array 16³"
    println("  Array 16³ vs flat_16: ✓  max=$(maximum(a16))  sum=$(sum(a16))")
end

if RUN_MARRAY
    println("  Specialising MArrays (may take a while — go get coffee)...")
    let an = zeros(UInt16, n, n, n), mn = fresh_mn()
        kw_array_n!(an, r_bench);  kw_marray_n!(mn, r_bench)
        @assert vec(an) == vec(mn) "marray_n index mapping differs from Array n³"
        println("  Array n³ vs marray_n:  ✓")
    end
    let a16 = zeros(UInt16, P, P, P), m16 = fresh_m16()
        kw_array_16!(a16, r_bench);  kw_marray_16!(m16, r_bench)
        @assert vec(a16) == vec(m16) "marray_16 index mapping differs from Array 16³"
        println("  Array 16³ vs marray_16: ✓")
    end
end
println("All sanity checks passed.\n")

# ─── Persistent read arrays (written once; used for all read benchmarks) ───────
const read_an  = let a=zeros(UInt16,n,n,n);   kw_array_n!(a,r_bench);  a end
const read_a16 = let a=zeros(UInt16,P,P,P);   kw_array_16!(a,r_bench); a end
const read_fn  = let f=zeros(UInt16,n^3);     kw_flat_n!(f,r_bench);   f end
const read_f16 = let f=zeros(UInt16,P^3);     kw_flat_16!(f,r_bench);  f end
if RUN_MARRAY
    const read_mn  = let m=fresh_mn();  kw_marray_n!(m,r_bench);  m end
    const read_m16 = let m=fresh_m16(); kw_marray_16!(m,r_bench); m end
end

# ─── WRITE benchmarks ─────────────────────────────────────────────────────────
# setup= reinitialises the target array per sample; not included in timing.
# This prevents UInt16 overflow across samples (each sample starts at zero).
println("══════════════════════════════════════════════════════")
println("WRITE: $(1000*998) read-modify-write ops per call")
println("══════════════════════════════════════════════════════\n")

println("=== 1w. Array{UInt16,3} $(n)³  (runtime stride=$n) ===")
display(@benchmark kw_array_n!(a, $r_bench) setup=(a=zeros(UInt16,n,n,n)))

println("\n=== 2w. Array{UInt16,3} $(P)³  (runtime stride=$P, power-of-2) ===")
display(@benchmark kw_array_16!(a, $r_bench) setup=(a=zeros(UInt16,P,P,P)))

println("\n=== 3w. flat Vector{UInt16} $(n^3), compile-time stride=$n ===")
display(@benchmark kw_flat_n!(f, $r_bench) setup=(f=zeros(UInt16,n^3)))

println("\n=== 4w. flat Vector{UInt16} $(P^3), bit-shift index (no multiply) ===")
display(@benchmark kw_flat_16!(f, $r_bench) setup=(f=zeros(UInt16,P^3)))

if RUN_MARRAY
    println("\n=== 5w. MArray{$(n)³, UInt16}  (StaticArrays compile-time strides) ===")
    display(@benchmark kw_marray_n!(m, $r_bench) setup=(m=fresh_mn()))

    println("\n=== 6w. MArray{$(P)³, UInt16}  (StaticArrays, power-of-2) ===")
    display(@benchmark kw_marray_16!(m, $r_bench) setup=(m=fresh_m16()))
end

# ─── READ benchmarks ──────────────────────────────────────────────────────────
println("\n\n══════════════════════════════════════════════════════")
println("READ: $(1000*998) reads per call  (primary backtracking pattern)")
println("══════════════════════════════════════════════════════\n")

println("=== 1r. Array{UInt16,3} $(n)³  (runtime stride=$n) ===")
display(@benchmark kr_array_n($read_an, $r_bench))

println("\n=== 2r. Array{UInt16,3} $(P)³  (runtime stride=$P, power-of-2) ===")
display(@benchmark kr_array_16($read_a16, $r_bench))

println("\n=== 3r. flat Vector{UInt16} $(n^3), compile-time stride=$n ===")
display(@benchmark kr_flat_n($read_fn, $r_bench))

println("\n=== 4r. flat Vector{UInt16} $(P^3), bit-shift index ===")
display(@benchmark kr_flat_16($read_f16, $r_bench))

if RUN_MARRAY
    println("\n=== 5r. MArray{$(n)³, UInt16}  (StaticArrays) ===")
    display(@benchmark kr_marray_n($read_mn, $r_bench))

    println("\n=== 6r. MArray{$(P)³, UInt16}  (StaticArrays, power-of-2) ===")
    display(@benchmark kr_marray_16($read_m16, $r_bench))
end

# ─── Max element (as originally requested) ────────────────────────────────────
println("\n─── Max element after one call ───")
let a=zeros(UInt16,n,n,n); kw_array_n!(a,r_bench)
    println("  Array $(n)³:       max = $(maximum(a))")
end
let f=zeros(UInt16,P^3); kw_flat_16!(f,r_bench)
    println("  flat Vector $(P)³: max = $(maximum(f))")
end

# ─── Native code inspection ────────────────────────────────────────────────────
# Uncomment to verify compiler choices. Key instructions:
#   imul / mulq      → runtime multiply (expected in cases 1, 2)
#   lea / shl / add  → strength reduction (expected in case 3, and case 5 if
#                      StaticArrays matches case 3)
#   shl / or         → pure bit-ops, no multiply at all (expected in cases 4, 6)
#
 let a=zeros(UInt16,n,n,n), f=zeros(UInt16,n^3), g=zeros(UInt16,P^3)
    println("\n--- 1. Array n³  (runtime multiply?) ---")
     @code_native kw_array_n!(a, r_bench)
     println("\n--- 3. flat n³   (strength reduction?) ---")
     @code_native kw_flat_n!(f, r_bench)
     println("\n--- 4. flat 16³  (shl/or only?) ---")
     @code_native kw_flat_16!(g, r_bench)
 end
 if RUN_MARRAY
     println("\n--- 5. MArray n³  (should match case 3) ---")
     @code_native kw_marray_n!(fresh_mn(), r_bench)
     println("\n--- 6. MArray 16³ (should match case 4) ---")
     @code_native kw_marray_16!(fresh_m16(), r_bench)
 end
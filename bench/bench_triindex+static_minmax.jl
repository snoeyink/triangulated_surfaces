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
using InteractiveUtils

const RUN_MARRAY = true   # ← set false to skip cases 5 & 6

const n = 13
const P = 16

Random.seed!(42)
const r_bench = rand(1:n, 1000)   # values in 1..13; fits all six array shapes

# ─── Helper functions for MArray setup ────────────────────────────────────────
# Defined before kernels so @benchmark setup= can call them.
@inline fresh_mn()  = (m = MArray{Tuple{n,n,n}, Int16}(undef); fill!(m, Int16(0)); m)
@inline fresh_m16() = (m = MArray{Tuple{P,P,P}, Int16}(undef); fill!(m, Int16(0)); m)

# ─── Write kernels ─────────────────────────────────────────────────────────────

# 1: Array n×n×n  — stride loaded from array header; LLVM uses runtime multiply
function kw_array_n!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += Int16(1)
    end
end

function kw_array_minmax_n!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        tri[mn, mx, r[i+2]] += Int16(1)
    end
end

# 2: Array 16×16×16 — same codegen path as case 1; different runtime stride
function kw_array_16!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        tri[r[i], r[i+1], r[i+2]] += Int16(1)
    end
end

function kw_array_minmax_16!(tri::Array{Int16,3}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        tri[mn, mx, r[i+2]] += Int16(1)
    end
end

# 3: flat Vector, stride n=13 is a compile-time constant.
#    LLVM strength-reduces: 13x = (8+4+1)x → shifts+adds, no imul.
function kw_flat_n!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) + n*(r[i+1]-1) + n*n*(r[i+2]-1)
        tri[idx+1] += Int16(1)
    end
end

function kw_flat_minmax_n!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        idx = (mn-1) + n*(mx-1) + n*n*(r[i+2]-1)
        tri[idx+1] += Int16(1)
    end
end

# 4: flat Vector 16³, pure bit-shift index — zero multiplies.
#    Valid because (x-1) ∈ 0..15 for x ∈ 1..16, so each field fits in 4 bits.
#    Equivalent to column-major Array{16,16,16} indexing (verified in sanity check).
function kw_flat_16!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) | ((r[i+1]-1) << 4) | ((r[i+2]-1) << 8)
        tri[idx+1] += Int16(1)
    end
end

function kw_flat_minmax_16!(tri::Vector{Int16}, r::Vector{Int})
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        idx = (mn-1) | ((mx-1) << 4) | ((r[i+2]-1) << 8)
        tri[idx+1] += Int16(1)
    end
end

# 5: MArray n×n×n — strides encoded in type parameters; LLVM sees constants.
#    Heap-allocated at this size, but no runtime stride lookup from header.
if RUN_MARRAY
    function kw_marray_n!(tri::MArray{Tuple{n,n,n}, Int16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            tri[r[i], r[i+1], r[i+2]] += Int16(1)
        end
    end

    function kw_marray_minmax_n!(tri::MArray{Tuple{n,n,n}, Int16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            mn,mx = minmax(r[i], r[i+1])
            tri[mn, mx, r[i+2]] += Int16(1)
        end
    end

    # 6: MArray 16×16×16 — power-of-2 strides, compile-time
    function kw_marray_16!(tri::MArray{Tuple{P,P,P}, Int16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            tri[r[i], r[i+1], r[i+2]] += Int16(1)
        end
    end

    function kw_marray_minmax_16!(tri::MArray{Tuple{P,P,P}, Int16}, r::Vector{Int})
        @inbounds for _ in 1:1000, i in 1:998
            mn,mx = minmax(r[i], r[i+1])
            tri[mn, mx, r[i+2]] += Int16(1)
        end
    end
end

# ─── Read kernels ──────────────────────────────────────────────────────────────
# Returns accumulated sum (prevents dead-code elimination; matches backtracking
# use pattern where tri_table is written once and read millions of times).

function kr_array_n(tri::Array{Int16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        s += tri[r[i], r[i+1], r[i+2]]
    end
    return s
end

function kr_array_minmax_n(tri::Array{Int16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        s += tri[mn, mx, r[i+2]]
    end
    return s
end

function kr_array_16(tri::Array{Int16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        s += tri[r[i], r[i+1], r[i+2]]
    end
    return s
end

function kr_array_minmax_16(tri::Array{Int16,3}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        s += tri[mn, mx, r[i+2]]
    end
    return s
end

function kr_flat_n(tri::Vector{Int16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) + n*(r[i+1]-1) + n*n*(r[i+2]-1)
        s += tri[idx+1]
    end
    return s
end

function kr_flat_minmax_n(tri::Vector{Int16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        idx = (mn-1) + n*(mx-1) + n*n*(r[i+2]-1)
        s += tri[idx+1]
    end
    return s
end

function kr_flat_16(tri::Vector{Int16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        idx = (r[i]-1) | ((r[i+1]-1) << 4) | ((r[i+2]-1) << 8)
        s += tri[idx+1]
    end
    return s
end

function kr_flat_minmax_16(tri::Vector{Int16}, r::Vector{Int})
    s = zero(Int32)
    @inbounds for _ in 1:1000, i in 1:998
        mn,mx = minmax(r[i], r[i+1])
        idx = (mn-1) | ((mx-1) << 4) | ((r[i+2]-1) << 8)
        s += tri[idx+1]
    end
    return s
end

if RUN_MARRAY
    function kr_marray_n(tri::MArray{Tuple{n,n,n}, Int16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            s += tri[r[i], r[i+1], r[i+2]]
        end
        return s
    end

    function kr_marray_minmax_n(tri::MArray{Tuple{n,n,n}, Int16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            mn,mx = minmax(r[i], r[i+1])
            s += tri[mn, mx, r[i+2]]
        end
        return s
    end

    function kr_marray_16(tri::MArray{Tuple{P,P,P}, Int16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            s += tri[r[i], r[i+1], r[i+2]]
        end
        return s
    end

    function kr_marray_minmax_16(tri::MArray{Tuple{P,P,P}, Int16}, r::Vector{Int})
        s = zero(Int32)
        @inbounds for _ in 1:1000, i in 1:998
            mn,mx = minmax(r[i], r[i+1])
            s += tri[mn, mx, r[i+2]]
        end
        return s
    end
end

# ─── Sanity checks ─────────────────────────────────────────────────────────────
println("Running sanity checks...")

# n-wide variants must agree on element values
let an = zeros(Int16, n, n, n),
    fn = zeros(Int16, n^3)
    kw_array_n!(an, r_bench);  kw_flat_n!(fn, r_bench)
    @assert vec(an) == vec(fn) "flat_n index mapping differs from Array n³"
    println("  Array n³ vs flat_n:  ✓  max=$(maximum(an))  sum=$(sum(an))")
end

let anm = zeros(Int16, n, n, n),
    fnm = zeros(Int16, n^3)
    kw_array_minmax_n!(anm, r_bench);  kw_flat_minmax_n!(fnm, r_bench)
    @assert vec(anm) == vec(fnm) "flat_minmax_n index mapping differs from Array n³ minmax"
    println("  Array n³ minmax vs flat_n minmax:  ✓  max=$(maximum(anm))  sum=$(sum(anm))")
end

# 16-wide variants must agree on element values
let a16 = zeros(Int16, P, P, P),
    f16 = zeros(Int16, P^3)
    kw_array_16!(a16, r_bench);  kw_flat_16!(f16, r_bench)
    @assert vec(a16) == vec(f16) "flat_16 index mapping differs from Array 16³"
    println("  Array 16³ vs flat_16: ✓  max=$(maximum(a16))  sum=$(sum(a16))")
end

let a16m = zeros(Int16, P, P, P),
    f16m = zeros(Int16, P^3)
    kw_array_minmax_16!(a16m, r_bench);  kw_flat_minmax_16!(f16m, r_bench)
    @assert vec(a16m) == vec(f16m) "flat_minmax_16 index mapping differs from Array 16³ minmax"
    println("  Array 16³ minmax vs flat_16 minmax: ✓  max=$(maximum(a16m))  sum=$(sum(a16m))")
end

if RUN_MARRAY
    println("  Specialising MArrays (may take a while — go get coffee)...")
    let an = zeros(Int16, n, n, n), mn = fresh_mn()
        kw_array_n!(an, r_bench);  kw_marray_n!(mn, r_bench)
        @assert vec(an) == vec(mn) "marray_n index mapping differs from Array n³"
        println("  Array n³ vs marray_n:  ✓")
    end
    let anm = zeros(Int16, n, n, n), mnm = fresh_mn()
        kw_array_minmax_n!(anm, r_bench);  kw_marray_minmax_n!(mnm, r_bench)
        @assert vec(anm) == vec(mnm) "marray_minmax_n index mapping differs from Array n³ minmax"
        println("  Array n³ minmax vs marray_n minmax:  ✓")
    end
    let a16 = zeros(Int16, P, P, P), m16 = fresh_m16()
        kw_array_16!(a16, r_bench);  kw_marray_16!(m16, r_bench)
        @assert vec(a16) == vec(m16) "marray_16 index mapping differs from Array 16³"
        println("  Array 16³ vs marray_16: ✓")
    end
    let a16m = zeros(Int16, P, P, P), m16m = fresh_m16()
        kw_array_minmax_16!(a16m, r_bench);  kw_marray_minmax_16!(m16m, r_bench)
        @assert vec(a16m) == vec(m16m) "marray_minmax_16 index mapping differs from Array 16³ minmax"
        println("  Array 16³ minmax vs marray_16 minmax: ✓")
    end
end
println("All sanity checks passed.\n")

# ─── Persistent read arrays (written once; used for all read benchmarks) ───────
const read_an  = let a=zeros(Int16,n,n,n);   kw_array_n!(a,r_bench);  a end
const read_an_mm = let a=zeros(Int16,n,n,n); kw_array_minmax_n!(a,r_bench); a end
const read_a16 = let a=zeros(Int16,P,P,P);   kw_array_16!(a,r_bench); a end
const read_a16_mm = let a=zeros(Int16,P,P,P); kw_array_minmax_16!(a,r_bench); a end
const read_fn  = let f=zeros(Int16,n^3);     kw_flat_n!(f,r_bench);   f end
const read_fn_mm = let f=zeros(Int16,n^3); kw_flat_minmax_n!(f,r_bench); f end
const read_f16 = let f=zeros(Int16,P^3);     kw_flat_16!(f,r_bench);  f end
const read_f16_mm = let f=zeros(Int16,P^3); kw_flat_minmax_16!(f,r_bench); f end
if RUN_MARRAY
    const read_mn  = let m=fresh_mn();  kw_marray_n!(m,r_bench);  m end
    const read_mn_mm = let m=fresh_mn(); kw_marray_minmax_n!(m,r_bench); m end
    const read_m16 = let m=fresh_m16(); kw_marray_16!(m,r_bench); m end
    const read_m16_mm = let m=fresh_m16(); kw_marray_minmax_16!(m,r_bench); m end
end

# ─── WRITE benchmarks ─────────────────────────────────────────────────────────
# setup= reinitialises the target array per sample; not included in timing.
# This prevents Int16 overflow across samples (each sample starts at zero).
println("══════════════════════════════════════════════════════")
println("WRITE: $(1000*998) read-modify-write ops per call")
println("══════════════════════════════════════════════════════\n")

println("=== 1w. Array{Int16,3} $(n)³  (runtime stride=$n) ===")
display(@benchmark kw_array_n!(a, $r_bench) setup=(a=zeros(Int16,n,n,n)))
println("=== 1w-mm. Array{Int16,3} $(n)³  (minmax on first two indices) ===")
display(@benchmark kw_array_minmax_n!(a, $r_bench) setup=(a=zeros(Int16,n,n,n)))

println("\n=== 2w. Array{Int16,3} $(P)³  (runtime stride=$P, power-of-2) ===")
display(@benchmark kw_array_16!(a, $r_bench) setup=(a=zeros(Int16,P,P,P)))
println("=== 2w-mm. Array{Int16,3} $(P)³  (minmax on first two indices) ===")
display(@benchmark kw_array_minmax_16!(a, $r_bench) setup=(a=zeros(Int16,P,P,P)))

println("\n=== 3w. flat Vector{Int16} $(n^3), compile-time stride=$n ===")
display(@benchmark kw_flat_n!(f, $r_bench) setup=(f=zeros(Int16,n^3)))
println("=== 3w-mm. flat Vector{Int16} $(n^3), minmax + compile-time stride=$n ===")
display(@benchmark kw_flat_minmax_n!(f, $r_bench) setup=(f=zeros(Int16,n^3)))

println("\n=== 4w. flat Vector{Int16} $(P^3), bit-shift index (no multiply) ===")
display(@benchmark kw_flat_16!(f, $r_bench) setup=(f=zeros(Int16,P^3)))
println("=== 4w-mm. flat Vector{Int16} $(P^3), minmax + bit-shift index ===")
display(@benchmark kw_flat_minmax_16!(f, $r_bench) setup=(f=zeros(Int16,P^3)))

if RUN_MARRAY
    println("\n=== 5w. MArray{$(n)³, Int16}  (StaticArrays compile-time strides) ===")
    display(@benchmark kw_marray_n!(m, $r_bench) setup=(m=fresh_mn()))
    println("=== 5w-mm. MArray{$(n)³, Int16}  (minmax on first two indices) ===")
    display(@benchmark kw_marray_minmax_n!(m, $r_bench) setup=(m=fresh_mn()))

    println("\n=== 6w. MArray{$(P)³, Int16}  (StaticArrays, power-of-2) ===")
    display(@benchmark kw_marray_16!(m, $r_bench) setup=(m=fresh_m16()))
    println("=== 6w-mm. MArray{$(P)³, Int16}  (minmax on first two indices) ===")
    display(@benchmark kw_marray_minmax_16!(m, $r_bench) setup=(m=fresh_m16()))
end

# ─── READ benchmarks ──────────────────────────────────────────────────────────
println("\n\n══════════════════════════════════════════════════════")
println("READ: $(1000*998) reads per call  (primary backtracking pattern)")
println("══════════════════════════════════════════════════════\n")

println("=== 1r. Array{Int16,3} $(n)³  (runtime stride=$n) ===")
display(@benchmark kr_array_n($read_an, $r_bench))
println("=== 1r-mm. Array{Int16,3} $(n)³  (minmax on first two indices) ===")
display(@benchmark kr_array_minmax_n($read_an_mm, $r_bench))

println("\n=== 2r. Array{Int16,3} $(P)³  (runtime stride=$P, power-of-2) ===")
display(@benchmark kr_array_16($read_a16, $r_bench))
println("=== 2r-mm. Array{Int16,3} $(P)³  (minmax on first two indices) ===")
display(@benchmark kr_array_minmax_16($read_a16_mm, $r_bench))

println("\n=== 3r. flat Vector{Int16} $(n^3), compile-time stride=$n ===")
display(@benchmark kr_flat_n($read_fn, $r_bench))
println("=== 3r-mm. flat Vector{Int16} $(n^3), minmax + compile-time stride=$n ===")
display(@benchmark kr_flat_minmax_n($read_fn_mm, $r_bench))

println("\n=== 4r. flat Vector{Int16} $(P^3), bit-shift index ===")
display(@benchmark kr_flat_16($read_f16, $r_bench))
println("=== 4r-mm. flat Vector{Int16} $(P^3), minmax + bit-shift index ===")
display(@benchmark kr_flat_minmax_16($read_f16_mm, $r_bench))

if RUN_MARRAY
    println("\n=== 5r. MArray{$(n)³, Int16}  (StaticArrays) ===")
    display(@benchmark kr_marray_n($read_mn, $r_bench))
    println("=== 5r-mm. MArray{$(n)³, Int16}  (minmax on first two indices) ===")
    display(@benchmark kr_marray_minmax_n($read_mn_mm, $r_bench))

    println("\n=== 6r. MArray{$(P)³, Int16}  (StaticArrays, power-of-2) ===")
    display(@benchmark kr_marray_16($read_m16, $r_bench))
    println("=== 6r-mm. MArray{$(P)³, Int16}  (minmax on first two indices) ===")
    display(@benchmark kr_marray_minmax_16($read_m16_mm, $r_bench))
end

# ─── Max element (as originally requested) ────────────────────────────────────
println("\n─── Max element after one call ───")
let a=zeros(Int16,n,n,n); kw_array_n!(a,r_bench)
    println("  Array $(n)³:       max = $(maximum(a))")
end
let f=zeros(Int16,P^3); kw_flat_16!(f,r_bench)
    println("  flat Vector $(P)³: max = $(maximum(f))")
end

# ─── Concise penalty summary (minmax / baseline) ─────────────────────────────
println("\n─── Penalty summary (minmax divided by baseline) ───")

t_w_array_n = @belapsed kw_array_n!(a, $r_bench) setup=(a=zeros(Int16,n,n,n))
t_w_array_n_mm = @belapsed kw_array_minmax_n!(a, $r_bench) setup=(a=zeros(Int16,n,n,n))
println("PENALTY write array_n ratio=$(round(t_w_array_n_mm / t_w_array_n, digits=3))")

t_w_array_16 = @belapsed kw_array_16!(a, $r_bench) setup=(a=zeros(Int16,P,P,P))
t_w_array_16_mm = @belapsed kw_array_minmax_16!(a, $r_bench) setup=(a=zeros(Int16,P,P,P))
println("PENALTY write array_16 ratio=$(round(t_w_array_16_mm / t_w_array_16, digits=3))")

t_w_flat_n = @belapsed kw_flat_n!(f, $r_bench) setup=(f=zeros(Int16,n^3))
t_w_flat_n_mm = @belapsed kw_flat_minmax_n!(f, $r_bench) setup=(f=zeros(Int16,n^3))
println("PENALTY write flat_n ratio=$(round(t_w_flat_n_mm / t_w_flat_n, digits=3))")

t_w_flat_16 = @belapsed kw_flat_16!(f, $r_bench) setup=(f=zeros(Int16,P^3))
t_w_flat_16_mm = @belapsed kw_flat_minmax_16!(f, $r_bench) setup=(f=zeros(Int16,P^3))
println("PENALTY write flat_16 ratio=$(round(t_w_flat_16_mm / t_w_flat_16, digits=3))")

t_r_array_n = @belapsed kr_array_n($read_an, $r_bench)
t_r_array_n_mm = @belapsed kr_array_minmax_n($read_an_mm, $r_bench)
println("PENALTY read array_n ratio=$(round(t_r_array_n_mm / t_r_array_n, digits=3))")

t_r_array_16 = @belapsed kr_array_16($read_a16, $r_bench)
t_r_array_16_mm = @belapsed kr_array_minmax_16($read_a16_mm, $r_bench)
println("PENALTY read array_16 ratio=$(round(t_r_array_16_mm / t_r_array_16, digits=3))")

t_r_flat_n = @belapsed kr_flat_n($read_fn, $r_bench)
t_r_flat_n_mm = @belapsed kr_flat_minmax_n($read_fn_mm, $r_bench)
println("PENALTY read flat_n ratio=$(round(t_r_flat_n_mm / t_r_flat_n, digits=3))")

t_r_flat_16 = @belapsed kr_flat_16($read_f16, $r_bench)
t_r_flat_16_mm = @belapsed kr_flat_minmax_16($read_f16_mm, $r_bench)
println("PENALTY read flat_16 ratio=$(round(t_r_flat_16_mm / t_r_flat_16, digits=3))")

if RUN_MARRAY
    t_w_marray_n = @belapsed kw_marray_n!(m, $r_bench) setup=(m=fresh_mn())
    t_w_marray_n_mm = @belapsed kw_marray_minmax_n!(m, $r_bench) setup=(m=fresh_mn())
    println("PENALTY write marray_n ratio=$(round(t_w_marray_n_mm / t_w_marray_n, digits=3))")

    t_w_marray_16 = @belapsed kw_marray_16!(m, $r_bench) setup=(m=fresh_m16())
    t_w_marray_16_mm = @belapsed kw_marray_minmax_16!(m, $r_bench) setup=(m=fresh_m16())
    println("PENALTY write marray_16 ratio=$(round(t_w_marray_16_mm / t_w_marray_16, digits=3))")

    t_r_marray_n = @belapsed kr_marray_n($read_mn, $r_bench)
    t_r_marray_n_mm = @belapsed kr_marray_minmax_n($read_mn_mm, $r_bench)
    println("PENALTY read marray_n ratio=$(round(t_r_marray_n_mm / t_r_marray_n, digits=3))")

    t_r_marray_16 = @belapsed kr_marray_16($read_m16, $r_bench)
    t_r_marray_16_mm = @belapsed kr_marray_minmax_16($read_m16_mm, $r_bench)
    println("PENALTY read marray_16 ratio=$(round(t_r_marray_16_mm / t_r_marray_16, digits=3))")
end

# ─── Native code inspection ────────────────────────────────────────────────────
# Uncomment to verify compiler choices. Key instructions:
#   imul / mulq      → runtime multiply (expected in cases 1, 2)
#   lea / shl / add  → strength reduction (expected in case 3, and case 5 if
#                      StaticArrays matches case 3)
#   shl / or         → pure bit-ops, no multiply at all (expected in cases 4, 6)
#
 let a=zeros(Int16,n,n,n), f=zeros(Int16,n^3), g=zeros(Int16,P^3)
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
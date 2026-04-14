using BenchmarkTools, Random

# ============================================================
# PARAMETERS
# ============================================================

const N = 13
const M = 1000

Random.seed!(42)

# independent random streams (harder for predictor)
const r1 = rand(1:N, M)
const r2 = rand(1:N, M)

# optional: permutation-based access (even less structure)
const perm = randperm(M)

# ============================================================
# PAIR INDEX VARIANTS
# ============================================================

@inline function pair_index_branch(i::Int, j::Int)::Int
    if i < j
        return (j-1)*(j-2) ÷ 2 + i
    else
        return (i-1)*(i-2) ÷ 2 + j
    end
end

@inline function pair_index_minmax(i::Int, j::Int)::Int
    mn, mx = minmax(i,j)
    return (mx-1)*(mx-2) ÷ 2 + mn
end

@inline function pair_index_ifelse(i::Int, j::Int)::Int
    cond = i < j
    a = ifelse(cond, i, j)
    b = ifelse(cond, j, i)
    return (b-1)*(b-2) ÷ 2 + a
end

@inline function pair_index_bit(i::Int, j::Int)::Int
    mask = -(i < j)
    a = (i & mask) | (j & ~mask)
    b = (j & mask) | (i & ~mask)
    return (b-1)*(b-2) ÷ 2 + a
end

const PAIR_LUT = let lut = Matrix{Int}(undef, N, N)
    for i in 1:N, j in 1:N
        if i == j
            lut[i,j] = 0
        elseif i < j
            lut[i,j] = (j-1)*(j-2) ÷ 2 + i
        else
            lut[i,j] = (i-1)*(i-2) ÷ 2 + j
        end
    end
    lut
end

@inline pair_index_lut(i::Int, j::Int)::Int = PAIR_LUT[i,j]

# ============================================================
# SANITY CHECK
# ============================================================

println("Running sanity checks...")
for i in 1:N, j in 1:N
    if i != j
        p = pair_index_branch(i,j)
        @assert p == pair_index_minmax(i,j)
        @assert p == pair_index_ifelse(i,j)
        @assert p == pair_index_bit(i,j)
        @assert p == pair_index_lut(i,j)
    end
end
println("All variants agree.\n")

# ============================================================
# BENCHMARK KERNELS
# ============================================================

# version A: independent random streams
function k_branch(r1, r2)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        s += pair_index_branch(r1[i], r2[i])
    end
    return s
end

function k_minmax(r1, r2)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        s += pair_index_minmax(r1[i], r2[i])
    end
    return s
end

function k_ifelse(r1, r2)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        s += pair_index_ifelse(r1[i], r2[i])
    end
    return s
end

function k_bit(r1, r2)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        s += pair_index_bit(r1[i], r2[i])
    end
    return s
end

function k_lut(r1, r2)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        s += pair_index_lut(r1[i], r2[i])
    end
    return s
end

# version B: permuted access (extra chaos)
function k_branch_perm(r1, r2, perm)
    s = 0
    @inbounds for _ in 1:1000, i in 1:M
        idx = perm[i]
        s += pair_index_branch(r1[idx], r2[idx])
    end
    return s
end

# ============================================================
# BENCHMARKS
# ============================================================

println("=== independent streams ===\n")

println("branch:")
display(@benchmark k_branch($r1, $r2))

println("\nminmax:")
display(@benchmark k_minmax($r1, $r2))

println("\nifelse:")
display(@benchmark k_ifelse($r1, $r2))

println("\nbit trick:")
display(@benchmark k_bit($r1, $r2))

println("\nlookup table:")
display(@benchmark k_lut($r1, $r2))

println("\n=== permuted access (branch stress test) ===\n")

println("branch (permuted):")
display(@benchmark k_branch_perm($r1, $r2, $perm))
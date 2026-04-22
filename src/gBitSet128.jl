module BitSets128

export BitSet128, BitSetOriented128, Rev
export singleton, valid_mask

# ======================================================================
# BitSet128
# ======================================================================

struct BitSet128
    words::NTuple{2, UInt64}
end

# Default empty constructor
BitSet128() = BitSet128((zero(UInt64), zero(UInt64)))

# Internal 0-based indexing (0:127)
@inline function singleton0(::Type{BitSet128}, i::Int)
    w1 = ifelse(UInt(i) < 64, one(UInt64) << i, zero(UInt64))
    w2 = ifelse(64 <= UInt(i) < 128, one(UInt64) << (i & 63), zero(UInt64))
    return BitSet128((w1, w2))
end

# Public 1-based indexing (1:128)
@inline singleton(::Type{BitSet128}, i::Int) = singleton0(BitSet128, i - 1)

@inline function valid_mask(::Type{BitSet128}, n::Int)
    n <= 0 && return BitSet128()
    n >= 128 && return BitSet128((typemax(UInt64), typemax(UInt64)))
    w1 = ifelse(n >= 64, typemax(UInt64), (one(UInt64) << n) - 1)
    w2 = ifelse(n <= 64, zero(UInt64), (one(UInt64) << (n - 64)) - 1)
    return BitSet128((w1, w2))
end

# ======================================================================
# BitSetOriented128
# ======================================================================

struct BitSetOriented128
    fwd::BitSet128
    rev::BitSet128
end

BitSetOriented128() = BitSetOriented128(BitSet128(), BitSet128())

@inline valid_mask(::Type{BitSetOriented128}, n::Int) = 
    BitSetOriented128(valid_mask(BitSet128, n), valid_mask(BitSet128, n))

# ======================================================================
# Base Overloads (BitSet128 & BitSetOriented128)
# ======================================================================

# Unary
Base.:~(a::BitSet128) = BitSet128((~a.words[1], ~a.words[2]))
Base.:~(a::BitSetOriented128) = BitSetOriented128(~a.fwd, ~a.rev)

# Binary bitwise
for op in (:|, :&, :xor)
    @eval Base.$op(a::BitSet128, b::BitSet128) = 
        BitSet128((Base.$op(a.words[1], b.words[1]), Base.$op(a.words[2], b.words[2])))
    
    @eval Base.$op(a::BitSetOriented128, b::BitSetOriented128) = 
        BitSetOriented128(Base.$op(a.fwd, b.fwd), Base.$op(a.rev, b.rev))
end

Base.setdiff(a::BitSet128, b::BitSet128) = a & (~b)
Base.setdiff(a::BitSetOriented128, b::BitSetOriented128) = BitSetOriented128(setdiff(a.fwd, b.fwd), setdiff(a.rev, b.rev))

# Comparisons & Queries
Base.:(==)(a::BitSet128, b::BitSet128) = a.words === b.words
Base.:(==)(a::BitSetOriented128, b::BitSetOriented128) = a.fwd == b.fwd && a.rev == b.rev

Base.isempty(a::BitSet128) = (a.words[1] | a.words[2]) == 0
Base.isempty(a::BitSetOriented128) = isempty(a.fwd) && isempty(a.rev)

Base.isdisjoint(a::BitSet128, b::BitSet128) = isempty(a & b)
Base.isdisjoint(a::BitSetOriented128, b::BitSetOriented128) = isdisjoint(a.fwd, b.fwd) && isdisjoint(a.rev, b.rev)

Base.in(i::Int, a::BitSet128) = !isdisjoint(a, singleton(BitSet128, i))
function Base.in(i::Int, a::BitSetOriented128)
    return i <= 128 ? in(i, a.fwd) : in(i - 128, a.rev)
end

# ======================================================================
# Iteration & Min Semantics
# ======================================================================

function Base.minimum(a::BitSet128)
    a.words[1] != 0 && return trailing_zeros(a.words[1]) + 1
    a.words[2] != 0 && return trailing_zeros(a.words[2]) + 65
    return 0
end

function Base.minimum(a::BitSetOriented128)
    m = minimum(a.fwd)
    m > 0 && return m
    m = minimum(a.rev)
    m > 0 && return m + 128
    return 0
end

# Iteration clears the lowest set bit at each step (BLSR instruction)
function Base.iterate(a::BitSet128, state = a.words)
    w1, w2 = state
    if w1 != 0
        return (trailing_zeros(w1) + 1, (w1 & (w1 - 1), w2))
    elseif w2 != 0
        return (trailing_zeros(w2) + 65, (zero(UInt64), w2 & (w2 - 1)))
    end
    return nothing
end

function Base.iterate(a::BitSetOriented128, state = (a.fwd.words..., a.rev.words...))
    fw1, fw2, rw1, rw2 = state
    if fw1 != 0
        return (trailing_zeros(fw1) + 1, (fw1 & (fw1 - 1), fw2, rw1, rw2))
    elseif fw2 != 0
        return (trailing_zeros(fw2) + 65, (zero(UInt64), fw2 & (fw2 - 1), rw1, rw2))
    elseif rw1 != 0
        return (trailing_zeros(rw1) + 129, (zero(UInt64), zero(UInt64), rw1 & (rw1 - 1), rw2))
    elseif rw2 != 0
        return (trailing_zeros(rw2) + 193, (zero(UInt64), zero(UInt64), zero(UInt64), rw2 & (rw2 - 1)))
    end
    return nothing
end

# ======================================================================
# Rev Wrapper
# ======================================================================

struct Rev
    x::BitSetOriented128
end

# Helper to swap wrapped inner elements
@inline _swap(x::BitSetOriented128) = BitSetOriented128(x.rev, x.fwd)

Base.:~(r::Rev) = ~_swap(r.x)

# Inject binary operations to cleanly intercept `Rev` semantics 
# where `Rev` triggers a pre-swap of its `fwd` and `rev` fields.
for op in (:|, :&, :xor, :setdiff, :(==), :isdisjoint)
    @eval Base.$op(a::BitSetOriented128, b::Rev) = Base.$op(a, _swap(b.x))
    @eval Base.$op(a::Rev, b::BitSetOriented128) = Base.$op(_swap(a.x), b)
    @eval Base.$op(a::Rev, b::Rev) = Base.$op(_swap(a.x), _swap(b.x))
end

end # module
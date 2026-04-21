module BitSetsOriented128

using Base: @propagate_inbounds

"""
    BitSetOriented128

Immutable set of directed edges on n ≤ 16 vertices (1-indexed), stored as two
`UInt128` words that share the same combinatorial bit-indexing scheme.

**Bit position** for undirected edge {i,j} with i < j (1-based):

    edge_pos(i,j) = i + (j-1)(j-2) ÷ 2      (max 120, for {15,16})

- `lo` bit at position p set  ⟺  forward edge i→j present
- `hi` bit at position p set  ⟺  reverse edge j→i present

Because both directions share the same bit position, reversal is free:

    reversed(s) = BitSetOriented128(s.hi, s.lo)

Use the lazy `Rev` wrapper for zero-overhead "op with reversed second argument":

    a ⊻ rev(b)   # XORs a with the reverse of b — no extra struct, no new op name
    a &  rev(b)
    a |  rev(b)
"""
struct BitSetOriented128
    lo::UInt128   # forward edges i→j  (i < j)
    hi::UInt128   # reverse edges j→i  (i < j), same bit positions as lo
end

BitSetOriented128() = BitSetOriented128(zero(UInt128), zero(UInt128))
Base.zero(::Type{BitSetOriented128}) = BitSetOriented128()

# ── Lazy reverse marker ─────────────────────────────────────────────────────────
# rev(b) as a second argument swaps b.lo ↔ b.hi at the call site, so the compiler
# sees plain UInt128 arithmetic — no heap allocation, no new operation name needed.

struct Rev
    s::BitSetOriented128
end

@inline rev(s::BitSetOriented128)      = Rev(s)          # lazy: O(0)
@inline reversed(s::BitSetOriented128) = BitSetOriented128(s.hi, s.lo)  # eager: O(1)

# ── Combinatorial coordinate ────────────────────────────────────────────────────

"""
    edge_pos(i, j) -> Int

1-based bit position for undirected edge {i,j}.  Caller must ensure `1 ≤ i < j ≤ 16`.
"""
@inline function edge_pos(i::Integer, j::Integer)::UInt8
    @boundscheck (1 ≤ i < j ≤ 16) ||
        throw(ArgumentError("vertices must be in increasing order from 1..16 (got $i, $j)"))
    UInt8(i) + (UInt8(j) - 0x01) * (UInt8(j) - 0x02) ÷ 0x02
end

# ── Singleton directed edge ─────────────────────────────────────────────────────

"""
    singleton(a, b) -> BitSetOriented128

Set containing only the directed edge a → b.
"""
@propagate_inbounds function singleton(a::Integer, b::Integer)
    i, j = minmax(a, b)
    @boundscheck (1 ≤ i < j ≤ 16) ||
        throw(ArgumentError("vertices must be distinct in 1..16 (got $a, $b)"))
    mask = one(UInt128) << (edge_pos(i, j) - 1)
    a < b ? BitSetOriented128(mask, zero(UInt128)) :
            BitSetOriented128(zero(UInt128), mask)
end

# ── Triangle ────────────────────────────────────────────────────────────────────

"""
    triangle(a, b, c) -> BitSetOriented128

Cyclic oriented triangle a→b, b→c, c→a.
Opposite winding: `reversed(triangle(a, b, c))`.
"""
@propagate_inbounds triangle(a::Integer, b::Integer, c::Integer) =
    singleton(a, b) | singleton(b, c) | singleton(c, a)

# ── Bitwise / set operations ────────────────────────────────────────────────────

@inline Base.:&(a::BitSetOriented128, b::BitSetOriented128) =
    BitSetOriented128(a.lo & b.lo, a.hi & b.hi)
@inline Base.:|(a::BitSetOriented128, b::BitSetOriented128) =
    BitSetOriented128(a.lo | b.lo, a.hi | b.hi)
@inline Base.:⊻(a::BitSetOriented128, b::BitSetOriented128) =
    BitSetOriented128(a.lo ⊻ b.lo, a.hi ⊻ b.hi)
@inline Base.:~(a::BitSetOriented128) =
    BitSetOriented128(~a.lo, ~a.hi)
    # Flips unused bits (positions ((n choose 2)+1)–128) too.
    # Safe in issubset / & contexts because the other operand keeps those bits 0.
    # For | and ⊻, the caller can use valid_mask(n) to reset extra bits to 0.

@inline function valid_mask(n::Integer)
    @boundscheck (0 ≤ n ≤ 16) ||
        throw(ArgumentError("n must be in 0..16, got $n"))
    k    = n * (n - 1) ÷ 2
    mask = (one(UInt128) << k) - one(UInt128)
    BitSetOriented128(mask, mask)
end
@inline complete_digraph(n::Integer) = valid_mask(n)

# ── Operations with lazy-reversed right-hand argument ───────────────────────────

@inline Base.:&(a::BitSetOriented128, b::Rev) =
    BitSetOriented128(a.lo & b.s.hi, a.hi & b.s.lo)
@inline Base.:|(a::BitSetOriented128, b::Rev) =
    BitSetOriented128(a.lo | b.s.hi, a.hi | b.s.lo)
@inline Base.:⊻(a::BitSetOriented128, b::Rev) =
    BitSetOriented128(a.lo ⊻ b.s.hi, a.hi ⊻ b.s.lo)
@inline Base.:~(a::Rev) = Rev(~a.s)

# ── Predicates ──────────────────────────────────────────────────────────────────

@inline Base.iszero(s::BitSetOriented128)  = iszero(s.lo) & iszero(s.hi)
@inline Base.isempty(s::BitSetOriented128) = iszero(s)
@inline Base.length(s::BitSetOriented128)  = count_ones(s.lo) + count_ones(s.hi)

@inline Base.issubset(a::BitSetOriented128, b::BitSetOriented128)   = iszero(a & ~b)
@inline Base.isdisjoint(a::BitSetOriented128, b::BitSetOriented128) = iszero(a & b)

# ── Minimum position ────────────────────────────────────────────────────────────

"""
    min_pos(s) -> UInt16

Index of the lowest set bit, with the encoding:
  •   1 – 128 → forward edge from `lo`, bit index = result
  • 129 – 256 → reverse edge from `hi`, bit index = result − 128

Returns 0 for an empty set.
"""
@inline function min_pos(s::BitSetOriented128)::UInt16
    iszero(s.lo) || return UInt16(trailing_zeros(s.lo)) + UInt16(1)
    iszero(s.hi) || return UInt16(trailing_zeros(s.hi)) + UInt16(129)
    UInt16(0)
end

# ── Conflict marking ─────────────────────────────────────────────────────────────

"""
    set_both(s, pos) -> BitSetOriented128

Set the bit at combinatorial position `pos` in **both** `lo` and `hi`, marking
both directions of that edge as present (e.g., a conflict after intersection).
"""
@inline function set_both(s::BitSetOriented128, pos::Integer)
    mask = one(UInt128) << (Int(pos) - 1)
    BitSetOriented128(s.lo | mask, s.hi | mask)
end

"""
    set_both(s, a, b) -> BitSetOriented128

Convenience form: mark both directions of edge {a, b}.
"""
@propagate_inbounds function set_both(s::BitSetOriented128, a::Integer, b::Integer)
    i, j = minmax(Int(a), Int(b))
    @boundscheck (1 ≤ i && j ≤ 16) ||
        throw(ArgumentError("vertices must be in 1..16"))
    set_both(s, edge_pos(i, j))
end

# ── Membership, add, remove ─────────────────────────────────────────────────────

"""
    in_directed(a, b, s) -> Bool

True if directed edge a→b is in `s`.
"""
@propagate_inbounds in_directed(a::Integer, b::Integer, s::BitSetOriented128) =
    !iszero(s & singleton(a, b))

"""
    push(s, a, b) -> BitSetOriented128

Return `s` with directed edge a→b added (non-mutating).
"""
@propagate_inbounds push(s::BitSetOriented128, a::Integer, b::Integer) =
    s | singleton(a, b)

"""
    delete(s, a, b) -> BitSetOriented128

Return `s` with directed edge a→b removed.
"""
@propagate_inbounds delete(s::BitSetOriented128, a::Integer, b::Integer) =
    s & ~singleton(a, b)

# ── Iteration ───────────────────────────────────────────────────────────────────
# Yields UInt16 positions in min_pos encoding: 1..128 (lo), 129..256 (hi).
# State is Tuple{UInt128,UInt128} — isbits, so no heap allocation.

Base.IteratorSize(::Type{BitSetOriented128}) = Base.HasLength()
Base.eltype(::Type{BitSetOriented128}) = UInt16

@inline function Base.iterate(s::BitSetOriented128,
                               state::Tuple{UInt128,UInt128} = (s.lo, s.hi))
    lo, hi = state
    if !iszero(lo)
        pos  = UInt16(trailing_zeros(lo)) + UInt16(1)
        rest = lo & (lo - one(UInt128))
        return (pos, (rest, hi))
    elseif !iszero(hi)
        pos  = UInt16(trailing_zeros(hi)) + UInt16(129)
        rest = hi & (hi - one(UInt128))
        return (pos, (lo, rest))
    else
        return nothing
    end
end

# ── Exports ──────────────────────────────────────────────────────────────────────

export BitSetOriented128, Rev,
       rev, reversed, edge_pos,
       singleton, triangle,
       set_both, push, delete, in_directed, min_pos, complete_digraph, valid_mask

end # module BitSetsOriented128

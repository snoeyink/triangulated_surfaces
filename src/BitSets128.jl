
#  from claude sonnet 4.6

export BitSet128, BitSetOriented128, Rev, singleton, valid_mask

# ─── BitSet128 ────────────────────────────────────────────────────────────────

"""
    BitSet128

Immutable bitset over public elements `1:128`.
`words[1]` holds internal bits 0–63 (elements 1–64);
`words[2]` holds internal bits 64–127 (elements 65–128).
All operations are allocation-free.
"""
struct BitSet128
    words::NTuple{2,UInt64}
end

BitSet128() = BitSet128((zero(UInt64), zero(UInt64)))
Base.eltype(::Type{BitSet128})    = Int
Base.IteratorSize(::Type{BitSet128}) = Base.SizeUnknown()

# ── internal 0-based helper (not exported) ────────────────────────────────────

@inline function singleton0(::Type{BitSet128}, i::Int)
    mask = one(UInt64) << (i & 63)
    i < 64 ? BitSet128((mask, zero(UInt64))) : BitSet128((zero(UInt64), mask))
end

# ── public constructors ───────────────────────────────────────────────────────

"""    singleton(BitSet128, i) → BitSet128 with only public element `i ∈ 1:128` set."""
@inline singleton(::Type{BitSet128}, i::Int) = singleton0(BitSet128, i - 1)

"""
    valid_mask(BitSet128, n) → BitSet128 with elements `1:n` set.

No masking is applied automatically elsewhere; call this explicitly when needed.
"""
function valid_mask(::Type{BitSet128}, n::Int)
    n <= 0   && return BitSet128()
    n >= 128 && return BitSet128((typemax(UInt64), typemax(UInt64)))
    if n <= 64
        w = n == 64 ? typemax(UInt64) : (one(UInt64) << n) - one(UInt64)
        return BitSet128((w, zero(UInt64)))
    else  # 65 ≤ n ≤ 127
        return BitSet128((typemax(UInt64), (one(UInt64) << (n - 64)) - one(UInt64)))
    end
end

# ── set operations ────────────────────────────────────────────────────────────

Base.:|(a::BitSet128, b::BitSet128)      = BitSet128(a.words .| b.words)
Base.:&(a::BitSet128, b::BitSet128)      = BitSet128(a.words .& b.words)
Base.xor(a::BitSet128, b::BitSet128)     = BitSet128(xor.(a.words, b.words))
Base.:~(a::BitSet128)                    = BitSet128(.~a.words)
Base.setdiff(a::BitSet128, b::BitSet128) = BitSet128(a.words .& .~b.words)

Base.:(==)(a::BitSet128, b::BitSet128)      = a.words == b.words
Base.isempty(a::BitSet128)                  = iszero(a.words[1] | a.words[2])
Base.isdisjoint(a::BitSet128, b::BitSet128) = isempty(a & b)
Base.issubset(a::BitSet128, b::BitSet128)   = isempty(setdiff(a, b))

function Base.in(i::Int, a::BitSet128)
    i0 = i - 1
    (i0 < 0 || i0 > 127) && return false
    @inbounds !iszero((a.words[(i0 >> 6) + 1] >> (i0 & 63)) & one(UInt64))
end

# ── minimum — returns 0 on empty (documented deviation from Base contract) ────

function Base.minimum(a::BitSet128)
    w1 = a.words[1];  !iszero(w1) && return trailing_zeros(w1) + 1
    w2 = a.words[2];  !iszero(w2) && return trailing_zeros(w2) + 65
    0
end

# ── iteration ─────────────────────────────────────────────────────────────────

@inline function _next128(w1::UInt64, w2::UInt64)
    if !iszero(w1)
        bit = trailing_zeros(w1)
        return (bit + 1,  (w1 & (w1 - one(UInt64)), w2))
    elseif !iszero(w2)
        bit = trailing_zeros(w2)
        return (bit + 65, (zero(UInt64), w2 & (w2 - one(UInt64))))
    else
        return nothing
    end
end

Base.iterate(a::BitSet128)                     = _next128(a.words[1], a.words[2])
Base.iterate(::BitSet128, s::NTuple{2,UInt64}) = _next128(s[1], s[2])


# ─── BitSetOriented128 ────────────────────────────────────────────────────────

"""
    BitSetOriented128

Immutable oriented bitset over public values `1:256`.
- Forward values `1:128` are stored in `fwd` (as `fwd` element `i`).
- Reverse values `129:256` are stored in `rev` (value `v` → `rev` element `v − 128`).

No automatic masking; iterators emit all set bits verbatim.
"""
struct BitSetOriented128
    fwd::BitSet128
    rev::BitSet128
end

BitSetOriented128() = BitSetOriented128(BitSet128(), BitSet128())
Base.eltype(::Type{BitSetOriented128})    = Int
Base.IteratorSize(::Type{BitSetOriented128}) = Base.SizeUnknown()

# ── public constructors ───────────────────────────────────────────────────────

"""    singleton(BitSetOriented128, i) → oriented singleton at value `i`."""
function singleton(::Type{BitSetOriented128}, i::Int)
    i <= 128 ?
        BitSetOriented128(singleton(BitSet128, i), BitSet128()) :
        BitSetOriented128(BitSet128(), singleton(BitSet128, i - 128))
end

"""    valid_mask(BitSetOriented128, n) → low-order `n` bits set in both `fwd` and `rev`."""
function valid_mask(::Type{BitSetOriented128}, n::Int)
    m = valid_mask(BitSet128, n)
    BitSetOriented128(m, m)
end

# ── set operations ────────────────────────────────────────────────────────────

Base.:|(a::BitSetOriented128, b::BitSetOriented128)      = BitSetOriented128(a.fwd | b.fwd, a.rev | b.rev)
Base.:&(a::BitSetOriented128, b::BitSetOriented128)      = BitSetOriented128(a.fwd & b.fwd, a.rev & b.rev)
Base.xor(a::BitSetOriented128, b::BitSetOriented128)     = BitSetOriented128(xor(a.fwd, b.fwd), xor(a.rev, b.rev))
Base.:~(a::BitSetOriented128)                            = BitSetOriented128(~a.fwd, ~a.rev)
Base.setdiff(a::BitSetOriented128, b::BitSetOriented128) = BitSetOriented128(setdiff(a.fwd, b.fwd), setdiff(a.rev, b.rev))

Base.:(==)(a::BitSetOriented128, b::BitSetOriented128)      = a.fwd == b.fwd && a.rev == b.rev
Base.isempty(a::BitSetOriented128)                          = isempty(a.fwd) && isempty(a.rev)
Base.isdisjoint(a::BitSetOriented128, b::BitSetOriented128) = isempty(a & b)
Base.issubset(a::BitSetOriented128, b::BitSetOriented128)   = isempty(setdiff(a, b))

function Base.in(i::Int, a::BitSetOriented128)
    1   <= i <= 128 && return (i       in a.fwd)
    129 <= i <= 256 && return ((i-128) in a.rev)
    false
end

# ── minimum — returns 0 on empty ──────────────────────────────────────────────

function Base.minimum(a::BitSetOriented128)
    m = minimum(a.fwd);  !iszero(m) && return m
    m = minimum(a.rev);  !iszero(m) && return m + 128
    0
end

# ── iteration ─────────────────────────────────────────────────────────────────
#
# Encoding of iterator offsets (all bit positions 0-based):
#   fwd.words[1] bits 0:63  →  oriented values   1: 64   (offset +1)
#   fwd.words[2] bits 0:63  →  oriented values  65:128   (offset +65)
#   rev.words[1] bits 0:63  →  oriented values 129:192   (offset +129)
#   rev.words[2] bits 0:63  →  oriented values 193:256   (offset +193)

@inline function _next_oriented(fw1::UInt64, fw2::UInt64, rw1::UInt64, rw2::UInt64)
    if !iszero(fw1)
        bit = trailing_zeros(fw1)
        return (bit + 1,   (fw1 & (fw1 - one(UInt64)), fw2, rw1, rw2))
    elseif !iszero(fw2)
        bit = trailing_zeros(fw2)
        return (bit + 65,  (zero(UInt64), fw2 & (fw2 - one(UInt64)), rw1, rw2))
    elseif !iszero(rw1)
        bit = trailing_zeros(rw1)
        return (bit + 129, (zero(UInt64), zero(UInt64), rw1 & (rw1 - one(UInt64)), rw2))
    elseif !iszero(rw2)
        bit = trailing_zeros(rw2)
        return (bit + 193, (zero(UInt64), zero(UInt64), zero(UInt64), rw2 & (rw2 - one(UInt64))))
    else
        return nothing
    end
end

function Base.iterate(a::BitSetOriented128)
    fw1, fw2 = a.fwd.words
    rw1, rw2 = a.rev.words
    _next_oriented(fw1, fw2, rw1, rw2)
end

Base.iterate(::BitSetOriented128, s::NTuple{4,UInt64}) =
    _next_oriented(s[1], s[2], s[3], s[4])


# ─── Rev ──────────────────────────────────────────────────────────────────────

"""
    Rev(x::BitSetOriented128)

Transparent wrapper presenting `x` with `fwd` and `rev` swapped when used as the
**final** operand in a binary operation:

    op(a, Rev(b))  ≡  op(a, BitSetOriented128(b.rev, b.fwd))

Supported unary: `~Rev(x)`.
Supported binary (Rev on right): `|`, `&`, `xor`, `setdiff`.
"""
struct Rev
    x::BitSetOriented128
end

@inline _unrev(r::Rev) = BitSetOriented128(r.x.rev, r.x.fwd)

Base.:~(r::Rev)                            = ~_unrev(r)
Base.:|(a::BitSetOriented128, r::Rev)      = a | _unrev(r)
Base.:&(a::BitSetOriented128, r::Rev)      = a & _unrev(r)
Base.xor(a::BitSetOriented128, r::Rev)     = xor(a, _unrev(r))
Base.setdiff(a::BitSetOriented128, r::Rev) = setdiff(a, _unrev(r))

#=


### Key design notes

| Concern            | Choice                                                                  |
| ------------------ | ----------------------------------------------------------------------- |
| Tuple update       | Ternary reconstruction in `singleton0`; `.op` broadcasting for bulk ops |
| Bit clearing       | `w & (w − 1)` — a single instruction, no separate mask                  |
| `minimum` on empty | Returns `0` by contract (documented deviation from `Base`)              |
| `Rev`              | Thin struct; `_unrev` inlines to zero overhead                          |
| No `length`        | `IteratorSize = SizeUnknown()` tells Julia not to expect it             |
| `@inbounds`        | Only where the index range is statically provable (`in`, `singleton0`)  |
| `singleton0`       | Not exported; `singleton` (1-based) is the public face                  |
| Oriented offsets   | `+1 / +65 / +129 / +193` exactly match the four word×half encoding      |

# Prompt

Use this as the new chat prompt:

Write a Julia module `BitSets128` using immutable wrapped tuples, concise and idiomatic.

Types:

* `struct BitSet128`

  * representation: `words::NTuple{2,UInt64}`
  * bits correspond internally to elements `0:127`
  * `words[1]` stores bits `0:63`, `words[2]` stores bits `64:127`
* `struct BitSetOriented128`

  * representation: `fwd::BitSet128`, `rev::BitSet128`
* also provide a wrapper `Rev` for `BitSetOriented128`

Public element conventions:

* `BitSet128` presents elements `1:128`
* internal helper `singleton0(BitSet128, i)` uses internal `0:127`
* public `singleton(BitSet128, i) = singleton0(BitSet128, i-1)`
* for “min-like” behavior, return `0` for empty rather than throwing
* `BitSetOriented128` presents values `1:255`

  * forward = `1:128`
  * reverse = `129:255`
  * iterator returns `Int`
  * iterator should return all nonzero bits; user is responsible for validity

Design/invariants:

* do not use `UInt128`
* do not use `StaticArrays`
* no mutation; use tuple reconstruction via `Base.setindex`
* no automatic masking of “invalid” oriented bits
* instead, for user masking when desired,  provide `valid_mask(BitSet128, n)` returning a`BitSet128` with the low-order `n` bits set (and `valid_mask(BitSetOriented128, n)` set both `fwd` and `rev`) 
* `BitSetOriented128` should allow all bits to exist; iterator/min operate on all nonzero bits

Operations to implement by overloading `Base` where appropriate:

* for `BitSet128` and `BitSetOriented128`:

  * `|`, `&`, `xor`, `~`, `setdiff`
  * `==`
  * `isempty`
  * `isdisjoint`
  * `in`
  * `iterate`
  * `minimum` or `min` semantics, but returning `0` on empty
* include the other standard BitSet-like conveniences that make sense, but do **not** provide `length`

Iteration/min semantics:

* `BitSet128`: public elements are `1:128`; internal bits are `0:127`
* `BitSetOriented128`: public values are `1:255`

  * iterate forward bits first, then reverse bits
  * no filtering of invalid bits; just emit all set bits according to the encoding

`Rev` wrapper:

* provide a wrapper `Rev(x::BitSetOriented128)`
* support `~Rev(x)` and binary ops involving `Rev` so that the last argument has `fwd` and `rev` swapped before the operation
* intent: convenient notation for “same operation, but with final argument reversed”
* implement this cleanly and idiomatically

Implementation preferences:

* concise, technical, idiomatic Julia
* wrapped-tuple value semantics
* use efficient bit primitives (`trailing_zeros`, `count_ones`, etc.) for iterator/min
* define internal helpers like `singleton0`
* export only the intended public API, not internal `0`-based helpers

Please include the full module code.
=#


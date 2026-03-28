# BitSet128.jl
# Compact bitset for elements 1..128 packed in a UInt128.
using Base: @propagate_inbounds

struct BitSet128
    bits::UInt128
end

BitSet128() = BitSet128(zero(UInt128))

@inline Base.:&(a::BitSet128, b::BitSet128)  = BitSet128(a.bits & b.bits)
@inline Base.:|(a::BitSet128, b::BitSet128)  = BitSet128(a.bits | b.bits)
@inline Base.:⊻(a::BitSet128, b::BitSet128)  = BitSet128(a.bits ⊻ b.bits)
@inline Base.:~(a::BitSet128)                = BitSet128(~a.bits)

@inline Base.iszero(a::BitSet128)   = iszero(a.bits)
@inline Base.length(a::BitSet128)   = count_ones(a.bits)

@inline function singleton(elem::Integer)
    @boundscheck (1 ≤ elem ≤ 128) || throw(ArgumentError("element $elem out of range [1, 128]"))
    BitSet128(UInt128(1) << (elem - 1))
end

@propagate_inbounds Base.in(elem::Integer, s::BitSet128)     = !iszero(s & singleton(elem))
@propagate_inbounds Base.push(s::BitSet128, elem::Integer)   = s | singleton(elem)
@propagate_inbounds Base.delete(s::BitSet128, elem::Integer) = s & ~singleton(elem)

@inline Base.issubset(a::BitSet128, b::BitSet128)   = iszero(a & ~b)
@inline Base.isdisjoint(a::BitSet128, b::BitSet128) = iszero(a & b)

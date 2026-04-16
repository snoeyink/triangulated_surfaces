module BitSets2n

using Base: @propagate_inbounds

_storage_type(::Val{8})   = UInt8
_storage_type(::Val{16})  = UInt16
_storage_type(::Val{32})  = UInt32
_storage_type(::Val{64})  = UInt64
_storage_type(::Val{128}) = UInt128

struct BitSet{N, T <: Unsigned}
    bits::T
end

@inline function BitSet{N}(bits::Integer) where {N}
    T = _storage_type(Val(N))
    BitSet{N, T}(T(bits))
end

BitSet{N}() where {N} = BitSet{N}(0)

const BitSet8   = BitSet{8,   UInt8}
const BitSet16  = BitSet{16,  UInt16}
const BitSet32  = BitSet{32,  UInt32}
const BitSet64  = BitSet{64,  UInt64}
const BitSet128 = BitSet{128, UInt128}

@inline Base.:&(a::BitSet{N,T}, b::BitSet{N,T}) where {N,T} = BitSet{N,T}(a.bits & b.bits)
@inline Base.:|(a::BitSet{N,T}, b::BitSet{N,T}) where {N,T} = BitSet{N,T}(a.bits | b.bits)
@inline Base.:⊻(a::BitSet{N,T}, b::BitSet{N,T}) where {N,T} = BitSet{N,T}(a.bits ⊻ b.bits)
@inline Base.:~(a::BitSet{N,T})                 where {N,T} = BitSet{N,T}(~a.bits)

@inline Base.iszero(a::BitSet) = iszero(a.bits)
@inline Base.isempty(a::BitSet) = iszero(a.bits)
@inline Base.length(a::BitSet) = count_ones(a.bits)

@propagate_inbounds @inline function singleton(::Type{BitSet{N,T}}, elem::Integer) where {N,T}
    @boundscheck (1 ≤ elem ≤ N) || throw(ArgumentError("element $elem out of range [1, $N]"))
    BitSet{N,T}(one(T) << (elem - 1))
end

@inline singleton(s::BitSet, elem::Integer) = singleton(typeof(s), elem)

@propagate_inbounds Base.in(elem::Integer, s::BitSet)     = !iszero(s & singleton(s, elem))
@propagate_inbounds push(s::BitSet, elem::Integer)        = s | singleton(s, elem)
@propagate_inbounds Base.delete(s::BitSet, elem::Integer) = s & ~singleton(s, elem)

@inline Base.issubset(a::BitSet{N,T},   b::BitSet{N,T}) where {N,T} = iszero(a & ~b)
@inline Base.isdisjoint(a::BitSet{N,T}, b::BitSet{N,T}) where {N,T} = iszero(a & b)

Base.IteratorSize(::Type{<:BitSet}) = Base.HasLength()

Base.eltype(::Type{<:BitSet{N}}) where {N} =
    N ≤ 255 ? UInt8 : UInt16 % insurance in case we ever add BitSet{256} or larger

@inline function Base.iterate(s::BitSet{N,T}, remaining::T = s.bits) where {N,T}
    iszero(remaining) && return nothing
    k    = eltype(s)(trailing_zeros(remaining)) + one(eltype(s))
    next = remaining & (remaining - one(T))
    return (k, next)
end

end # module BitSets2n
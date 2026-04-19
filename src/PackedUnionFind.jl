module PackedUnionFind

export PackedUF, find_root, union_sets

# 0xFEDCBA9876543210 assigns parent[i] = i for i in 0:15 in 4-bit chunks
const INIT_STATE = 0xFEDCBA9876543210

struct PackedUF
    state::UInt64
end

PackedUF() = PackedUF(INIT_STATE)

"""
    find_root(uf::PackedUF, i::Int)
Returns the 1-based root of the 1-based element `i`.
"""
@inline function find_root(uf::PackedUF, i::Int)
    curr = i - 1 # 0-based for bit shifting
    @inbounds while true
        # Read the 4-bit chunk corresponding to 'curr'
        p = (uf.state >> (curr << 2)) & 0x0F
        p == curr && return curr + 1
        curr = p
    end
end

"""
    union_sets(uf::PackedUF, i::Int, j::Int)
Attempts to union elements `i` and `j`.
Returns `(new_uf, true)` if successful, or `(uf, false)` if a cycle is detected.
"""
@inline function union_sets(uf::PackedUF, i::Int, j::Int)
    root_i = find_root(uf, i) - 1
    root_j = find_root(uf, j) - 1
    
    # Premature cycle detected in the star!
    root_i == root_j && return uf, false
    
    # Set parent of root_j to root_i
    mask = ~(UInt64(0x0F) << (root_j << 2))
    new_state = (uf.state & mask) | (UInt64(root_i) << (root_j << 2))
    
    return PackedUF(new_state), true
end

end # module
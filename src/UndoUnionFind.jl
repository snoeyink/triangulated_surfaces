module UndoUnionFind

export StarMatrix, find_root, union_sets!, undo!

struct StarMatrix
    # parents[neighbor, vertex] = parent of 'neighbor' in the link of 'vertex'
    parents::Matrix{Int8}
end

function StarMatrix(N::Int)
    parents = Matrix{Int8}(undef, N, N)
    for v in 1:N, n in 1:N
        parents[n, v] = n # Initialize parent[i] = i
    end
    return StarMatrix(parents)
end

@inline function find_root(sm::StarMatrix, vertex::Int, i::Int)
    curr = Int8(i)
    @inbounds while true
        p = sm.parents[curr, vertex]
        p == curr && return curr
        curr = p
    end
end

"""
    union_sets!(sm::StarMatrix, vertex::Int, i::Int, j::Int)
Unions `i` and `j` in the star of `vertex`. 
Returns `(success, undo_record)`.
"""
@inline function union_sets!(sm::StarMatrix, vertex::Int, i::Int, j::Int)
    root_i = find_root(sm, vertex, i)
    root_j = find_root(sm, vertex, j)
    
    # Cycle detected
    root_i == root_j && return false, (0, 0, Int8(0))
    
    # Record the old state for undo
    @inbounds old_parent = sm.parents[root_j, vertex]
    
    # Mutate
    @inbounds sm.parents[root_j, vertex] = root_i
    
    return true, (vertex, root_j, old_parent)
end

"""
    undo!(sm::StarMatrix, record::Tuple{Int, Int, Int8})
Restores the exact pointer modified by `union_sets!`.
"""
@inline function undo!(sm::StarMatrix, record::Tuple{Int, Int, Int8})
    vertex, child, old_parent = record
    if vertex != 0
        @inbounds sm.parents[child, vertex] = old_parent
    end
    return nothing
end

end # module
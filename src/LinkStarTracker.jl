module LinkStarTracker

export LinkTracker, add_triangle!, undo!

struct LinkTracker
    degree::Vector{UInt8}
    star_next::Vector{UInt8}
end

"""
    LinkTracker(N::UInt8, max_edge_idx::UInt8)

"""
function LinkTracker(N, max_edge_idx)
    # Both arrays fit easily in CPU cache (e.g., 16 + 210 bytes)
    return LinkTracker(zeros(UInt8, N), zeros(UInt8, max_edge_idx))
end

"""
    check_cycle(tracker, vertex, start_node, next_node, edge_index)

Traverses the link of `vertex` starting from `next_node`.
Returns:
   0 if it forms an open chain.
   1 if it forms a valid closed loop (size == degree).
  -1 if it forms a premature cycle (size < degree).
"""
@inline function check_cycle(tracker::LinkTracker, vertex, start_node, next_node, edge_index)
    curr = next_node
    len = 1
    @inbounds while true
        nxt = tracker.star_next[edge_index(UInt8(vertex), UInt8(curr))]
        nxt == 0 && return 0 # Reached the end of an open chain
        
        curr = nxt
        len += 1
        
        if curr == start_node
            # We joined two ends! Check if it includes all triangles at this vertex.
            return len == tracker.degree[vertex] ? 1 : -1
        end
    end
end

"""
    add_triangle!(tracker, i, j, k, edge_index)

Records the addition of triangle (i, j, k).
Returns `(success, closed_i, closed_j, closed_k)`.
If `success == false`, a premature cycle was detected, and you MUST still call `undo!`.
"""
@inline function add_triangle!(tracker::LinkTracker, i, j, k, edge_index)
    i, j, k = UInt8(i), UInt8(j), UInt8(k)
    @inbounds begin
        # 1. Update degrees
        tracker.degree[i] += 1
        tracker.degree[j] += 1
        tracker.degree[k] += 1
        
        # 2. Add the link edges (intrusive linked list pointers)
        tracker.star_next[edge_index(i, j)] = k
        tracker.star_next[edge_index(j, k)] = i
        tracker.star_next[edge_index(k, i)] = j
    end
    
    # 3. Traverse to see if we closed a loop equal to the degree
    c_i = check_cycle(tracker, i, j, k, edge_index)
    c_j = check_cycle(tracker, j, k, i, edge_index)
    c_k = check_cycle(tracker, k, i, j, edge_index)
    
    if c_i == -1 || c_j == -1 || c_k == -1
        return false, false, false, false
    end
    
    return true, c_i == 1, c_j == 1, c_k == 1
end

"""
    undo!(tracker, i, j, k, edge_index)

Perfect lazy deletion: restores the exact state before `add_triangle!` was called.
"""
@inline function undo!(tracker::LinkTracker, i, j, k, edge_index)
    i, j, k = UInt8(i), UInt8(j), UInt8(k)
    @inbounds begin
        tracker.degree[i] -= 1
        tracker.degree[j] -= 1
        tracker.degree[k] -= 1
        
        tracker.star_next[edge_index(i,j)] = 0
        tracker.star_next[edge_index(j,k)] = 0
        tracker.star_next[edge_index(k,i)] = 0
    end
    return nothing
end

end # module
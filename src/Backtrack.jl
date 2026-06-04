

function enumerate_triangulated_surfaces(points::Vector{Point3D})
    n = length(points)
    !(MIN_VERTICES <= n <= MAX_VERTICES) && throw(ArgumentError("Need $(MIN_VERTICES)<= n <= $(MAX_VERTICES) points"))
    n == N || throw(ArgumentError("This implementation currently requires n == N == $(N)"))

    _, emap, econflT = precompute_conflicts(points)
    nv_complete = UInt8(0)
    has = BitSetOriented128()
    confl = BitSet128()
    vlinks = [PackedUnionFind.PackedUF() for _ in 1:N]
    
    out = BitSetOriented128[]
    count = Ref{Int64}(0) 
    
    # Create plot
    pl = create_plot(points)



    function backtrack!(::Val{NV}, has::BitSetOriented128, confl::BitSet128,
                        out::Vector{BitSetOriented128}, count::Base.RefValue{Int64}) where {NV}
        
        @inline uconfl(e::UInt8) = (u_edge(e)) in confl

        @inline function updatevlink(v::UInt8, a::UInt8, b::UInt8)
            next_state, ok = PackedUnionFind.union_sets(vlinks[v], a, b)
            if !ok # a,b already in same component of the link of v. 
                PackedUnionFind.single_component(next_state, v) | return false # link of v cannot be completed to a manifold
                nv_complete += UInt8(1) # completed the link of v
                @inbounds confl |= usesV[v] 
            end
            @inbounds vlinks[v] = next_state # 
            return true
        end

        bdry = setdiff(has, Rev(has))
        if isempty(bdry)
            if nv_complete == NV-3 # need to count completed vertices 
                count[] += 1
                push!(out, has)
            end
            return
        end

        e = minimum(bdry)

        if D
            e == 0 && error("Can't have minimum edge e==0") 
            uconfl(e) && error("Can't have minimum edge e in conflict set")
        end
        
        @inbounds i, j = emap[e]
        for k::UInt8 in 1:NV
            (k == i || k == j) && continue

            eki = edge_index(k, i)
            uconfl(eki) && continue
            ejk = edge_index(j, k)
            uconfl(ejk) && continue
            tindex = triangle_index(i, j, k)
            @inbounds (isdisjoint(econflT[tindex], has.fwd) && isdisjoint(econflT[tindex], has.rev)) || continue
            
            # Add triangles
            plot_tri(pl, i, j, k)
            @inbounds begin           
                save_vi = vlinks[i]
                save_vj = vlinks[j]
                save_vk = vlinks[k]
                save_nv_complete = nv_complete
                save_confl = confl # this is because I'm using the conflict set to forbid edges to completed vertices, but I could do this another way.

                if updatevlink(k, j, i) && 
                    updatevlink(i, k, j) && 
                    updatevlink(j, i, k)

                    backtrack!(Val(NV), 
                        has | singleton(BitSetOriented128, eki) | singleton(BitSetOriented128, ejk), 
                        confl | econflT[tindex], out, count)
                end
                vlinks[i] = save_vi
                vlinks[j] = save_vj
                vlinks[k] = save_vk
                nv_complete = save_nv_complete
                confl = save_confl
                plot_pop!(pl)
            end
        end
    end

    @inline function seed_vlinks!(vlinks::Vector{PackedUnionFind.PackedUF}, i::UInt8, j::UInt8, k::UInt8, l::UInt8)
        @inbounds begin
            vlinks[i], ok = PackedUnionFind.union_sets(vlinks[i], j, k); ok || return false
            vlinks[i], ok = PackedUnionFind.union_sets(vlinks[i], j, l); ok || return false

            vlinks[j], ok = PackedUnionFind.union_sets(vlinks[j], i, k); ok || return false
            vlinks[j], ok = PackedUnionFind.union_sets(vlinks[j], i, l); ok || return false

            vlinks[k], ok = PackedUnionFind.union_sets(vlinks[k], i, j); ok || return false
            vlinks[l], ok = PackedUnionFind.union_sets(vlinks[l], i, j); ok || return false
        end
        return true
    end


    for i in 1:(N-1), j in (i+1):N
        maxe = edge_index(i, j)
        higher = setdiff(valid_mask(BitSet128, 128), valid_mask(BitSet128, maxe))

        for k in 1:N-1
            (k == i || k == j) && continue
            ejk = edge_index(j, k)
            u_edge(ejk) < maxe || continue
            eki = edge_index(k, i)
            u_edge(eki) < maxe || continue
            hasijk = singleton(BitSetOriented128, edge_index(i, j)) |
                    singleton(BitSetOriented128, edge_index(j, k)) |
                    singleton(BitSetOriented128, edge_index(k, i)) |
                    singleton(BitSetOriented128, edge_index(j, i)) 
            t1 = triangle_index(i, j, k)
            conflk = higher | econflT[t1]
            (isdisjoint(hasijk.fwd, conflk) && isdisjoint(hasijk.rev, conflk)) || continue

            for l in k+1:N
                (l == i || l == j) && continue
             
                t2 = triangle_index(j, i, l)
                has = hasijk |
                      singleton(BitSetOriented128, edge_index(i, l)) |
                      singleton(BitSetOriented128, edge_index(l, j))
                
                confl = conflk | econflT[t2]
                (isdisjoint(has.fwd, confl) && isdisjoint(has.rev, confl)) || continue

                seed_vlinks!(vlinks, UInt8(i), UInt8(j), UInt8(k), UInt8(l)) || continue

                backtrack!(Val(N), has, confl, out, count)
            end
        end
    end

    return count[], out
end
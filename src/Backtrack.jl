@inline reverse_oriented(s::BitSetOriented128) = BitSetOriented128(s.rev, s.fwd)

@inline function make_oriented_triangle(i::Int, j::Int, k::Int)
    tri = singleton(BitSetOriented128, edge_index(i, j)) |
          singleton(BitSetOriented128, edge_index(j, k)) |
          singleton(BitSetOriented128, edge_index(k, i))
    return reverse_oriented(tri)
end

@inline function updatevlink!(vlinks::Vector{PackedUnionFind.PackedUF}, v::Int, a::Int, b::Int)
    next_state, ok = PackedUnionFind.union_sets(vlinks[v], a, b)
    ok || return false
    @inbounds vlinks[v] = next_state
    return true
end

@inline is_complete(has::BitSetOriented128) = isempty(setdiff(has, Rev(has)))

function backtrack!(::Val{NV}, has::BitSetOriented128, confl::BitSet128,
                    vlinks::Vector{PackedUnionFind.PackedUF},
                    econflT::Vector{BitSet128}, emap::Vector{EdgeIJ},
                    out::Vector{BitSetOriented128}, count::Base.RefValue{Int}) where {NV}
    bdry = setdiff(has, Rev(has))
    if isempty(bdry)
        count[] += 1
        push!(out, has)
        return
    end

    e = minimum(bdry)
    e == 0 && return
    ue = undirected_edge_index(e)
    (ue in confl) && return

    @inbounds i, j = emap[e]
    for k in 1:NV
        (k == i || k == j) && continue

        eki = edge_index(k, i)
        ueki = undirected_edge_index(eki)
        (ueki in confl) && continue

        ejk = edge_index(j, k)
        uejk = undirected_edge_index(ejk)
        (uejk in confl) && continue

        @inbounds begin
            save_vi = vlinks[i]
            save_vj = vlinks[j]
            save_vk = vlinks[k]

            if updatevlink!(vlinks, k, j, i) &&
               updatevlink!(vlinks, i, k, j) &&
               updatevlink!(vlinks, j, i, k)
                tri = make_oriented_triangle(i, j, k)
                tindex = triangle_index(i, j, k)
                backtrack!(Val(NV), has | tri, confl | econflT[tindex], vlinks, econflT, emap, out, count)
            end

            vlinks[i] = save_vi
            vlinks[j] = save_vj
            vlinks[k] = save_vk
        end
    end
end

@inline function seed_vlinks!(vlinks::Vector{PackedUnionFind.PackedUF}, i::Int, j::Int, k::Int, l::Int)
    ok = true

    ok || return false
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

@inline function choose_seed_orientation(has::BitSetOriented128)
    rev_has = reverse_oriented(has)
    minimum(has) <= minimum(rev_has) ? has : rev_has
end

function enumerate_triangulated_surfaces(points::Vector{Point3D})
    n = length(points)
    !(MIN_VERTICES <= n <= MAX_VERTICES) && throw(ArgumentError("Need $(MIN_VERTICES)<= n <= $(MAX_VERTICES) points"))
    n == N || throw(ArgumentError("This implementation currently requires n == N == $(N)"))

    _, emap, econflT = precompute_conflicts(points)
    out = BitSetOriented128[]
    count = Ref(0)

    for i in 1:(N-1), j in (i+1):N
        maxf = edge_index(i, j)
        higher = setdiff(valid_mask(BitSet128, 128), valid_mask(BitSet128, maxf))

        for k in 1:N
            (k == i || k == j) && continue
            for l in 1:N
                (l == i || l == j || l == k) && continue

                t1 = triangle_index(i, j, k)
                t2 = triangle_index(j, i, l)

                has = singleton(BitSetOriented128, edge_index(i, j)) |
                      singleton(BitSetOriented128, edge_index(j, k)) |
                      singleton(BitSetOriented128, edge_index(k, i)) |
                      singleton(BitSetOriented128, edge_index(j, i)) |
                      singleton(BitSetOriented128, edge_index(i, l)) |
                      singleton(BitSetOriented128, edge_index(l, j))
                has = choose_seed_orientation(has)

                confl = higher | econflT[t1] | econflT[t2]

                vlinks = [PackedUnionFind.PackedUF() for _ in 1:N]
                seed_vlinks!(vlinks, i, j, k, l) || continue

                backtrack!(Val(N), has, confl, vlinks, econflT, emap, out, count)
            end
        end
    end

    return count[], out
end
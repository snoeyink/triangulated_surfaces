# ─────────────────────────────────────────────────────────────────────────────
# Backtrack.jl
# Requires: BdryLoop.jl, BitSet128.jl, EdgesTriangles.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── SearchFrame ───────────────────────────────────────────────────────────────
# isbits: stored inline in Vector{SearchFrame}, no boxing.
# Field order minimises padding: UInt128-aligned fields first.
# sizeof: 16+16+4+1+1+1+1+1 = 41 → padded to 48 by allocator (16-byte alignment)
struct SearchFrame
    added_edgeset    ::BitSet128   # b.added_edgeset before the move
    forbidden_edgeset::BitSet128   # b.forbidden_edgeset before the move
    ve               ::BdryVE      # b.v[vertex_i] before the move
    n_unused         ::Int8        # b.n_unused before the move
    vertex_i         ::Int8        # boundary vertex modified (also serves as debug check)
    resume_i         ::Int8        # 0 = level exhausted; else: resume from this vertex
    resume_k         ::Int8        # 0 = try link at resume_i; k>0 = try ears from k
    saved_head       ::Int8        # b.head before this move (restored by popBE!)
end

# ── BdryStack ─────────────────────────────────────────────────────────────────
mutable struct BdryStack{D}
    frames ::Vector{SearchFrame}
    sp     ::Int
end

function BdryStack(; debug::Bool = true)
    BdryStack{debug}(
        Vector{SearchFrame}(undef, MAX_VERTICES),
        0,
    )
end

# ── Conflict checks ───────────────────────────────────────────────────────────

# ear (i, j, k=nxt(i)): 2 new edges (i,j) and (j,k)
@inline function ear_ok(b                ::BdryLoop,
                         tri_table       ::Array{Int16,3},
                         edgesets::Vector{Tri_Edgesets},
                         i::Int, j::Int) :: Bool
    @inbounds begin
        ti = tri_table[i, j, b.v[i].next]
        return isdisjoint(edgesets[ti].conf, b.added_edgeset) && isdisjoint(edgesets[ti].has, b.forbidden_edgeset)
    end
end

# link (i, j=nxt(i), k=nxt2(i)): 1 new edge (i,k)
@inline function link_ok(b                ::BdryLoop,
                         tri_table        ::Array{Int16,3},
                         edgesets::Vector{Tri_Edgesets},
                         i::Int) :: Bool
    @inbounds begin
        j  = Int(b.v[i].next)
        ti = tri_table[i, j, b.v[j].next]
        return isdisjoint(edgesets[ti].conf, b.added_edgeset) && isdisjoint(edgesets[ti].has, b.forbidden_edgeset)
    end
end

# ── Min-removable check ───────────────────────────────────────────────────────
# Returns true iff b.v[i].tri ≤ every other removable boundary vertex's tri.
# Call AFTER the move that set b.v[i].tri.
@inline function is_min_removable(b::BdryLoop, i::Int) :: Bool
    t = @inbounds b.v[i].tri
    m = b.head
    @inbounds while true
        if removable(b, m) && b.v[m].tri < t
            return false
        end
        m = Int(b.v[m].next)
        m == b.head && break
    end
    return true
end

# ── Resume helpers ────────────────────────────────────────────────────────────
@inline function first_unused(b::BdryLoop) :: Int
    @inbounds for k in 1:b.n; b.status[k] == UNUSED && return k; end
    return 0
end

@inline function next_unused_after(b::BdryLoop, after::Int) :: Int
    @inbounds for k in after+1:b.n; b.status[k] == UNUSED && return k; end
    return 0
end

# After link(i): next sibling is ear(i, first_unused), or advance to nxt(i)
@inline function resume_after_link(b::BdryLoop, i::Int) :: Tuple{Int8,Int8}
    k = first_unused(b)
    k > 0 && return Int8(i), Int8(k)
    ni = nxt(b, i)
    ni == b.head && return Int8(0), Int8(0)   # wrapped: level exhausted
    return Int8(ni), Int8(0)
end

# After ear(i, k): next sibling is ear(i, k') or advance to nxt(i)
@inline function resume_after_ear(b::BdryLoop, i::Int, k::Int) :: Tuple{Int8,Int8}
    k2 = next_unused_after(b, k)
    k2 > 0 && return Int8(i), Int8(k2)
    ni = nxt(b, i)
    ni == b.head && return Int8(0), Int8(0)
    return Int8(ni), Int8(0)
end

# ── add_ear! ──────────────────────────────────────────────────────────────────
# Prereqs: on_boundary(b,i), unused(b,k); ear_ok() checked by caller.
# Pushes SearchFrame; resume_i/resume_k encode the next sibling at THIS level.
function add_ear!(b                ::BdryLoop,
                  s                ::BdryStack{D},
                  tri_table        ::Array{Int16,3},
                  edgesets::Vector{Tri_Edgesets},
                  i::Int, k::Int,
                  resume_i::Int8, resume_k::Int8) where {D}
    D && @assert on_boundary(b, i) "add_ear!: i not on boundary"
    D && @assert unused(b, k)      "add_ear!: k not unused"
    @inbounds begin
        j = Int(b.v[i].next)
        t = tri_table[i, j, k]
        D && @assert !iszero(Int(t)) "add_ear!: triangle not in tri_table"

        # Save full state before any modification
        s.sp += 1
        s.frames[s.sp] = SearchFrame(
            b.added_edgeset, b.forbidden_edgeset,
            b.v[i], Int8(b.n_unused), Int8(i),
            resume_i, resume_k, Int8(b.head),
        )

        # Topology
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)
        b.v[k]      = BdryVE(Int8(j), Int8(i), t)
        b.status[k] = ON_BOUNDARY
        b.n_unused -= 1

        # Edgesets: 2 new edges (j,k) and (i,k)
        b.added_edgeset     |= edgesets[t].has
        b.forbidden_edgeset |= edgesets[t].conf
    end
    return b
end

# ── add_link! ─────────────────────────────────────────────────────────────────
# Prereqs: on_boundary(b,i), on_boundary(nxt(b,i)); link_ok() checked by caller.
# v[j] preserved untouched (needed for undo). b.head updated if j was head.
function add_link!(b                ::BdryLoop,
                   s                ::BdryStack{D},
                   tri_table        ::Array{Int16,3},
                   edgesets::Vector{Tri_Edgesets},
                   i::Int,
                   resume_i::Int8, resume_k::Int8) where {D}
    D && @assert on_boundary(b, i)          "add_link!: i not on boundary"
    D && @assert on_boundary(b, nxt(b, i))  "add_link!: nxt(i) not on boundary"
    @inbounds begin
        j = Int(b.v[i].next)
        k = Int(b.v[j].next)
        t = tri_table[i, j, k]
        D && @assert !iszero(Int(t)) "add_link!: triangle not in tri_table"

        # Save full state before any modification
        s.sp += 1
        s.frames[s.sp] = SearchFrame(
            b.added_edgeset, b.forbidden_edgeset,
            b.v[i], Int8(b.n_unused), Int8(i),
            resume_i, resume_k, Int8(b.head),
        )

        # Topology: j becomes interior; v[j] preserved for undo
        b.status[j] = INTERIOR
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)

        # Keep head on the boundary
        if j == b.head; b.head = k; end

        # Edgesets: 1 new edge (i,k)
        b.added_edgeset     |= edgesets[t].has
        b.forbidden_edgeset |= edgesets[t].conf
    end
    return b
end

# ── popBE! ────────────────────────────────────────────────────────────────────
# Undo the last pushed move for vertex i.
# Restores topology, edgesets, n_unused, head from the frame.
function popBE!(b::BdryLoop, s::BdryStack{D}, i::Int) where {D}
    D && @assert s.sp > 0      "popBE!: stack empty"
    D && @assert !ear1(b, i)   "popBE!: precondition !ear1 failed"
    @inbounds begin
        frame = s.frames[s.sp]
        D && @assert Int(frame.vertex_i) == i "popBE!: wrong vertex"

        if ear0(b, i)
            k = Int(b.v[i].next)
            D && @assert on_boundary(b, k) "popBE!: ear0: k not on_boundary"
            b.status[k] = UNUSED
        else                               # link(i)
            j = Int(b.v[i].opp)
            D && @assert interior(b, j)    "popBE!: link: j not interior"
            b.status[j] = ON_BOUNDARY
        end

        b.v[i]              = frame.ve
        b.n_unused          = Int(frame.n_unused)
        b.added_edgeset     = frame.added_edgeset
        b.forbidden_edgeset = frame.forbidden_edgeset
        b.head              = Int(frame.saved_head)
        s.sp -= 1
    end
    return b
end


# ── Backtracking skeleton ─────────────────────────────────────────────────────
# Advance loop (per level):
#   for i around the boundary from (resume_i, resume_k):
#     if resume_k==0 && link(b,i) && removable(b,i) && link_ok(b,T,E,i):
#       ri,rk = resume_after_link(b,i)
#       add_link!(b,s,T,E,i,ri,rk)
#       is_min_removable(b,i) ? go deeper : popBE!(b,s,i)
#     for k = max(resume_k,1)..b.n:
#       if unused(b,k) && ear_ok(b,T,E,i,k):
#         ri,rk = resume_after_ear(b,i,k)
#         add_ear!(b,s,T,E,i,k,ri,rk)
#         is_min_removable(b,i) ? go deeper : popBE!(b,s,i)
#     resume_k = 0 (fresh for subsequent vertices)
# Terminal: b.n_unused==0 && nxt3(b,head)==head → check last triangle → write
# Backtrack: no advance → pop frame → restore → jump to (frame.resume_i, frame.resume_k)
#
# TODO: implement after data structure tests pass
function backtrack!(b  ::BdryLoop,
                    s  ::BdryStack{D},
                    T  ::Array{Int16,3},
                    es ::Vector{Tri_Edgesets},
                    out::IO) where {D}
    error("backtrack!: not yet implemented")
end
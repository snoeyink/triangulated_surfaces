# ==== Small isbits structs ====

struct BdryVE
    next::Int8
    opp::Int8
    tri::Int16
end

# ─────────────────────────────────────────────────────────────────────────────
# BoundaryLoop.jl
# External: triangle_index(i,j,k)::Int16  (i<j<k required)
#           t_index(i,j,k)::Int16          (sorts internally)
# ─────────────────────────────────────────────────────────────────────────────

# ── BdryVE ───────────────────────────────────────────────────────────────────
# sizeof == 4, isbits, no padding.
# Vertex i owns directed edge i → next.
# Pending left triangle is (i, next, opp) with index tri.
# Int8 safe: n ≤ 16, vertex indices 1..16 fit.
struct BdryVE
    next::Int8
    opp ::Int8
    tri ::Int16
end

# ── Vertex status ─────────────────────────────────────────────────────────────
const UNUSED      = Int8(0)
const ON_BOUNDARY = Int8(1)
const INTERIOR    = Int8(2)

# ── BdryLoop ──────────────────────────────────────────────────────────────────
# D::Bool is resolved at compile time:
#   D = true  → assertions + push-index bookkeeping active
#   D = false → all debug code compiled away; dbg_idx has length 0
# Use BdryLoop(n) for development, BdryLoop(n; debug=false) for production.
mutable struct BdryLoop{D}
    v       ::Vector{BdryVE}   # half-edge data, 1..n
    status  ::Vector{Int8}     # UNUSED / ON_BOUNDARY / INTERIOR
    stack   ::Vector{BdryVE}   # undo stack, pre-allocated length n
    dbg_idx ::Vector{Int8}     # [D=true] vertex index at each push level
    sp      ::Int              # stack pointer (0 = empty)
    n       ::Int              # vertex capacity
    head    ::Int              # a vertex currently in the boundary loop
end

function BdryLoop(n::Int; debug::Bool = true)
    BdryLoop{debug}(
        Vector{BdryVE}(undef, n),
        fill(UNUSED, n),
        Vector{BdryVE}(undef, n),
        debug ? Vector{Int8}(undef, n) : Int8[],
        0, n, 0,
    )
end

# ── Accessors (always @inbounds; return Int for arithmetic) ───────────────────
@inline nxt(b::BdryLoop, i::Int)  = @inbounds Int(b.v[i].next)
@inline opp(b::BdryLoop, i::Int)  = @inbounds Int(b.v[i].opp)
@inline tri(b::BdryLoop, i::Int)  = @inbounds b.v[i].tri
@inline nxt2(b::BdryLoop, i::Int) = nxt(b, nxt(b, i))
@inline nxt3(b::BdryLoop, i::Int) = nxt(b, nxt2(b, i))

# ── Status predicates ─────────────────────────────────────────────────────────
@inline unused(b::BdryLoop, i::Int)      = @inbounds b.status[i] == UNUSED
@inline on_boundary(b::BdryLoop, i::Int) = @inbounds b.status[i] == ON_BOUNDARY
@inline interior(b::BdryLoop, i::Int)    = @inbounds b.status[i] == INTERIOR

# ── Boundary vertex predicates ────────────────────────────────────────────────
@inline ear0(b::BdryLoop, i::Int) = opp(b, i) == nxt2(b, i)
@inline ear1(b::BdryLoop, i::Int) = nxt(b, opp(b, i)) == i
@inline link(b::BdryLoop, i::Int) = !ear0(b, i) && !ear1(b, i)

# Removable: safe to add the pending triangle (no non-manifold edge created).
# Not removable iff link(i) AND opp(i) is already in the boundary loop.
@inline removable(b::BdryLoop, i::Int) = !link(b, i) || !on_boundary(b, opp(b, i))

@inline function first_tri(b::BdryLoop, i::Int)
    j = nxt(b, i); k = nxt(b, j)
    ear0(b, i) && ear1(b, i) && i < j && i < k
end

@inline function is_last(b::BdryLoop, i::Int)
    nxt3(b, i) == i && i < nxt(b, i) && i < nxt2(b, i)
end

# ── init_loop!(b, i, j, k) ───────────────────────────────────────────────────
# Requires i < j < k. Resets all vertex statuses.
function init_loop!(b::BdryLoop{D}, i::Int, j::Int, k::Int) where {D}
    D && @assert i < j < k "init_loop!: need i<j<k"
    t = triangle_index(i, j, k)
    @inbounds begin
        fill!(b.status, UNUSED)
        b.v[i] = BdryVE(Int8(j), Int8(k), t)
        b.v[j] = BdryVE(Int8(k), Int8(i), t)
        b.v[k] = BdryVE(Int8(i), Int8(j), t)
        b.status[i] = ON_BOUNDARY
        b.status[j] = ON_BOUNDARY
        b.status[k] = ON_BOUNDARY
    end
    b.sp = 0; b.head = i
    return b
end

# ── add_ear!(b, i, k) ────────────────────────────────────────────────────────
# Glue ear triangle (i, nxt(i), k); unused vertex k joins the boundary.
# Pushes v[i]. Post: nxt(i)=k, opp(i)=j; ear0(i)=true, ear1(k)=true.
function add_ear!(b::BdryLoop{D}, i::Int, k::Int) where {D}
    D && @assert on_boundary(b, i) "add_ear!: i not on boundary"
    D && @assert unused(b, k)      "add_ear!: k not unused"
    @inbounds begin
        j = Int(b.v[i].next)
        t = t_index(i, j, k)
        b.sp += 1
        b.stack[b.sp] = b.v[i]
        if D; b.dbg_idx[b.sp] = Int8(i); end
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)
        b.v[k]      = BdryVE(Int8(j), Int8(i), t)
        b.status[k] = ON_BOUNDARY
    end
    return b
end

# ── add_link!(b, i) ──────────────────────────────────────────────────────────
# Glue link triangle (i, j, k) where j=nxt(i), k=nxt2(i).
# j becomes INTERIOR; v[j] is preserved untouched (implicit undo storage).
# Pushes v[i]. Post: nxt(i)=k, opp(i)=j (interior); link(i)=true.
#
#    Must be BdryVE(k,j,t): v[j].next=k, so BdryVE(j,k,t) would give
#    opp(i)=k and nxt(j)=k=opp → ear0(i)=true, killing the link post-condition
#    and sending popBE! into the wrong branch.
function add_link!(b::BdryLoop{D}, i::Int) where {D}
    D && @assert on_boundary(b, i)          "add_link!: i not on boundary"
    D && @assert on_boundary(b, nxt(b, i))  "add_link!: nxt(i) not on boundary"
    @inbounds begin
        j = Int(b.v[i].next)
        k = Int(b.v[j].next)
        t = t_index(i, j, k)
        b.sp += 1
        b.stack[b.sp] = b.v[i]
        if D; b.dbg_idx[b.sp] = Int8(i); end
        b.status[j] = INTERIOR           # v[j] preserved for undo
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)
    end
    return b
end

# ── popBE!(b, i) ─────────────────────────────────────────────────────────────
# Undo the operation that last pushed for vertex i.
# ear0 branch: nxt(i) was added by add_ear! → mark it UNUSED.
# link branch: opp(i) was buried by add_link! → restore it to ON_BOUNDARY.
# Restores v[i] from the stack.
function popBE!(b::BdryLoop{D}, i::Int) where {D}
    D && @assert b.sp > 0                  "popBE!: stack empty"
    D && @assert !ear1(b, i)               "popBE!: precondition !ear1 failed"
    D && @assert Int(b.dbg_idx[b.sp]) == i "popBE!: stack top is wrong vertex"
    @inbounds begin
        if ear0(b, i)
            k = Int(b.v[i].next)
            D && @assert on_boundary(b, k) "popBE!: ear0 branch: k not on_boundary"
            b.status[k] = UNUSED
        else                               # !ear0 && !ear1 ⟹ link(i)
            j = Int(b.v[i].opp)
            D && @assert interior(b, j)    "popBE!: link branch: j not interior"
            b.status[j] = ON_BOUNDARY
        end
        b.v[i] = b.stack[b.sp]
        b.sp  -= 1
    end
    return b
end


function removable(vc, a, nbdry)   
        bv = falses(16)
        for i = 1:nbdry
            bv[a] = true
            a = vc[a].next
        end
        mn = triangle_index(14,15,16)+1 # larger than any triangle index
        mi = a
        for i = 1:nbdry
            if vc[a].tri < mn && (!bv[vc[a].opp] || e1(a))
               mn = vc[a].tri
               mi = a
            end
            a = vc[a].next
        end
        return mn, mi
    end
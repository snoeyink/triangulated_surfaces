# ─────────────────────────────────────────────────────────────────────────────
# BdryLoop.jl
# External: 
#           t_index(i,j,k)::Int16  
# ─────────────────────────────────────────────────────────────────────────────
function t_index(a::Int, b::Int, c::Int)
    mn = min(a,b,c)
    mx = max(a,b,c)
    md = a-mn+c-mx+b
    return triangle_index(mn, md, mx)  # +1 for 1-based indexing
end

const MAX_VERTICES = 16

# ── BdryVE ───────────────────────────────────────────────────────────────────
# sizeof == 4, isbits, no padding.
# Vertex i owns directed edge i → next.
# Pending left triangle is (i, next, opp) with index tri.
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
# D=true  → assertions active (compiled away when D=false)
mutable struct BdryLoop{D}
    v      ::Vector{BdryVE}   # half-edge data, 1..MAX_VERTICES
    status ::Vector{Int8}     # UNUSED / ON_BOUNDARY / INTERIOR
    head   ::Int              # a vertex currently on the boundary
end

function BdryLoop(; debug::Bool = true)
    BdryLoop{debug}(
        Vector{BdryVE}(undef, MAX_VERTICES),
        fill(UNUSED, MAX_VERTICES),
        0,
    )
end

# ── BdryStack ─────────────────────────────────────────────────────────────────
# Separate undo stack for backtracking.
# D=true  → push-index bookkeeping active (compiled away when D=false)
mutable struct BdryStack{D}
    data    ::Vector{BdryVE}
    dbg_idx ::Vector{Int8}    # [D=true] vertex index at each push level
    sp      ::Int             # stack pointer (0 = empty)
end

function BdryStack(; debug::Bool = true)
    BdryStack{debug}(
        Vector{BdryVE}(undef, MAX_VERTICES),
        debug ? Vector{Int8}(undef, MAX_VERTICES) : Int8[],
        0,
    )
end

# ── Accessors ────────────────────────────────────────────────────────────────
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

@inline removable(b::BdryLoop, i::Int) = !link(b, i) || !on_boundary(b, opp(b, i))

@inline function first_tri(b::BdryLoop, i::Int)
    j = nxt(b, i); k = nxt(b, j)
    ear0(b, i) && ear1(b, i) && i < j && i < k
end

# ── init_loop!(b, i, j, k) ───────────────────────────────────────────────────
# Requires i < j < k. Resets all vertex statuses.
function init_loop!(b::BdryLoop{D}, i::Int, j::Int, k::Int) where {D}
    D && @assert i < j < k "init_loop!: need i<j<k"
    t = t_index(i, j, k)
    @inbounds begin
        fill!(b.status, UNUSED)
        b.v[i] = BdryVE(Int8(j), Int8(k), t)
        b.v[j] = BdryVE(Int8(k), Int8(i), t)
        b.v[k] = BdryVE(Int8(i), Int8(j), t)
        b.status[i] = ON_BOUNDARY
        b.status[j] = ON_BOUNDARY
        b.status[k] = ON_BOUNDARY
    end
    b.head = i
    return b
end

# ── add_ear!(b, s, i, k) ─────────────────────────────────────────────────────
# Glue ear triangle (i, nxt(i), k); unused vertex k joins the boundary.
# Pushes v[i] to s. Post: nxt(i)=k, opp(i)=j; ear0(i)=true, ear1(k)=true.
function add_ear!(b::BdryLoop{D}, s::BdryStack{D}, i::Int, k::Int) where {D}
    D && @assert on_boundary(b, i) "add_ear!: i not on boundary"
    D && @assert unused(b, k)      "add_ear!: k not unused"
    @inbounds begin
        j = Int(b.v[i].next)
        t = t_index(i, j, k)
        s.sp += 1
        s.data[s.sp] = b.v[i]
        if D; s.dbg_idx[s.sp] = Int8(i); end
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)
        b.v[k]      = BdryVE(Int8(j), Int8(i), t)
        b.status[k] = ON_BOUNDARY
    end
    return b
end

# ── add_link!(b, s, i) ───────────────────────────────────────────────────────
# Glue link triangle (i, j, k) where j=nxt(i), k=nxt2(i).
# j becomes INTERIOR; v[j] is preserved untouched (implicit undo storage).
# Pushes v[i] to s. Post: nxt(i)=k, opp(i)=j (interior); link(i)=true.
function add_link!(b::BdryLoop{D}, s::BdryStack{D}, i::Int) where {D}
    D && @assert on_boundary(b, i)         "add_link!: i not on boundary"
    D && @assert on_boundary(b, nxt(b, i)) "add_link!: nxt(i) not on boundary"
    @inbounds begin
        j = Int(b.v[i].next)
        k = Int(b.v[j].next)
        t = t_index(i, j, k)
        s.sp += 1
        s.data[s.sp] = b.v[i]
        if D; s.dbg_idx[s.sp] = Int8(i); end
        b.status[j] = INTERIOR            # v[j] preserved for undo
        b.v[i]      = BdryVE(Int8(k), Int8(j), t)
    end
    return b
end

# ── popBE!(b, s, i) ──────────────────────────────────────────────────────────
# Undo the last push for vertex i.
# ear0 branch: nxt(i) was added by add_ear! → mark it UNUSED.
# link branch: opp(i) was buried by add_link! → restore it to ON_BOUNDARY.
function popBE!(b::BdryLoop{D}, s::BdryStack{D}, i::Int) where {D}
    D && @assert s.sp > 0                  "popBE!: stack empty"
    D && @assert !ear1(b, i)               "popBE!: precondition !ear1 failed"
    D && @assert Int(s.dbg_idx[s.sp]) == i "popBE!: stack top is wrong vertex"
    @inbounds begin
        if ear0(b, i)
            k = Int(b.v[i].next)
            D && @assert on_boundary(b, k) "popBE!: ear0 branch: k not on_boundary"
            b.status[k] = UNUSED
        else                               # link(i)
            j = Int(b.v[i].opp)
            D && @assert interior(b, j)    "popBE!: link branch: j not interior"
            b.status[j] = ON_BOUNDARY
        end
        b.v[i] = s.data[s.sp]
        s.sp  -= 1
    end
    return b
end
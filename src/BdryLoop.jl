# ─────────────────────────────────────────────────────────────────────────────
# BdryLoop.jl
# External: triangle_index(i,j,k)::Int16  requires i<j<k
# ─────────────────────────────────────────────────────────────────────────────

const MAX_VERTICES = 16

struct BdryVE
    next::Int8
    opp ::Int8
    tri ::Int16
end

const UNUSED      = Int8(0)
const ON_BOUNDARY = Int8(1)
const INTERIOR    = Int8(2)

# D=true: assertions active; compiled away entirely when D=false
mutable struct BdryLoop{D}
    v                ::Vector{BdryVE}   # half-edge data, length MAX_VERTICES
    status           ::Vector{Int8}     # UNUSED / ON_BOUNDARY / INTERIOR
    head             ::Int              # any vertex currently on the boundary
    n                ::Int              # total vertices for this search instance (fixed)
    n_unused         ::Int              # current count of UNUSED vertices
    added_edgeset    ::BitSet128        # edges of every triangle added so far
    forbidden_edgeset::BitSet128        # edges conflicting with any added triangle
end

function BdryLoop(n::Int; debug::Bool = true)
    BdryLoop{debug}(
        Vector{BdryVE}(undef, MAX_VERTICES),
        fill(UNUSED, MAX_VERTICES),
        0, n, n,
        BitSet128(), BitSet128(),
    )
end

# ── Accessors ────────────────────────────────────────────────────────────────
@inline nxt(b::BdryLoop, i::Int)  = @inbounds Int(b.v[i].next)
@inline opp(b::BdryLoop, i::Int)  = @inbounds Int(b.v[i].opp)
@inline tri(b::BdryLoop, i::Int)  = @inbounds b.v[i].tri
@inline nxt2(b::BdryLoop, i::Int) = nxt(b, nxt(b, i))
@inline nxt3(b::BdryLoop, i::Int) = nxt(b, nxt2(b, i))

# ── Status predicates ────────────────────────────────────────────────────────
@inline unused(b::BdryLoop, i::Int)      = @inbounds b.status[i] == UNUSED
@inline on_boundary(b::BdryLoop, i::Int) = @inbounds b.status[i] == ON_BOUNDARY
@inline interior(b::BdryLoop, i::Int)    = @inbounds b.status[i] == INTERIOR

# ── Boundary predicates ──────────────────────────────────────────────────────
@inline ear0(b::BdryLoop, i::Int)      = opp(b, i) == nxt2(b, i)
@inline ear1(b::BdryLoop, i::Int)      = nxt(b, opp(b, i)) == i
@inline link(b::BdryLoop, i::Int)      = !ear0(b, i) && !ear1(b, i)
@inline removable(b::BdryLoop, i::Int) = !link(b, i) || !on_boundary(b, opp(b, i))

@inline function first_tri(b::BdryLoop, i::Int)
    j = nxt(b, i); k = nxt(b, j)
    ear0(b, i) && ear1(b, i) && i < j && i < opp(b, i)
end

# ── init_loop! (test version) ─────────────────────────────────────────────────
# Uses triangle_index; zeroes edgesets.  Tests depend on this form.
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
    b.head              = i
    b.n_unused          = b.n - 3
    b.added_edgeset     = BitSet128()
    b.forbidden_edgeset = BitSet128()
    return b
end

# ── init_loop! (search version) ───────────────────────────────────────────────
# Takes renumbered triangle index t and precomputed edgesets from the main loop.
function init_loop!(b::BdryLoop{D}, i::Int, j::Int, k::Int, t::Int16,
                    init_added::BitSet128, init_forbidden::BitSet128) where {D}
    D && @assert i < j < k "init_loop!: need i<j<k"
    @inbounds begin
        fill!(b.status, UNUSED)
        b.v[i] = BdryVE(Int8(j), Int8(k), t)
        b.v[j] = BdryVE(Int8(k), Int8(i), t)
        b.v[k] = BdryVE(Int8(i), Int8(j), t)
        b.status[i] = ON_BOUNDARY
        b.status[j] = ON_BOUNDARY
        b.status[k] = ON_BOUNDARY
    end
    b.head              = i
    b.n_unused          = b.n - 3
    b.added_edgeset     = init_added
    b.forbidden_edgeset = init_forbidden
    return b
end
# ─────────────────────────────────────────────────────────────────────────────
# BdryLoop.jl
# External dependencies (must be included/loaded before this file):
#   BitSet128.jl   — provides BitSet128
#   triangle_index(i,j,k)::Int16  requires i < j < k
# ─────────────────────────────────────────────────────────────────────────────

const MAX_VERTICES = 16

# ── BdryVE ───────────────────────────────────────────────────────────────────
# sizeof == 4, isbits, no padding.
# Vertex i owns directed boundary edge i → next.
# Pending left triangle is (i, next, opp) with renumbered index tri.
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
# D=true  → debug assertions active
# D=false → assertions compiled away entirely
mutable struct BdryLoop
    v                ::Vector{BdryVE}  # half-edge data, lenght n
    status           ::Vector{Int8}    # UNUSED / ON_BOUNDARY / INTERIOR
    head             ::Int             # any vertex currently on the boundary
    n                ::Int             # total vertices (fixed for this search)
    n_unused         ::Int             # current count of UNUSED vertices
    added_edgeset    ::BitSet128       # edges of every triangle added so far
    forbidden_edgeset::BitSet128       # edges conflicting with any added triangle
end

function BdryLoop(n::Int)
    BdryLoop(
        Vector{BdryVE}(undef, n),
        fill(UNUSED, n),
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

# ── Status predicates ─────────────────────────────────────────────────────────
@inline unused(b::BdryLoop, i::Int)      = @inbounds b.status[i] == UNUSED
@inline on_boundary(b::BdryLoop, i::Int) = @inbounds b.status[i] == ON_BOUNDARY
@inline interior(b::BdryLoop, i::Int)    = @inbounds b.status[i] == INTERIOR

# ── Boundary vertex predicates ────────────────────────────────────────────────
@inline ear0(b::BdryLoop, i::Int) = opp(b, i) == nxt2(b, i)
@inline ear1(b::BdryLoop, i::Int) = nxt(b, opp(b, i)) == i
@inline link(b::BdryLoop, i::Int) = !ear0(b, i) && !ear1(b, i)

# Removable: safe to add the pending triangle (no non-manifold edge).
# Not removable iff link(i) AND opp(i) is already on the boundary.
@inline removable(b::BdryLoop, i::Int) = !link(b, i) || !on_boundary(b, opp(b, i))

# first_tri: true for exactly one vertex when the boundary is a single triangle.
@inline function first_tri(b::BdryLoop, i::Int)
    j = nxt(b, i); k = nxt(b, j)
    ear0(b, i) && ear1(b, i) && i < j && i < opp(b, i)
end

# ── init_loop! (test version) ─────────────────────────────────────────────────
# Uses triangle_index directly; leaves edgesets zeroed.
# Call this form from tests and anywhere renumbering is not needed.
function init_loop!(b::BdryLoop, i::Int, j::Int, k::Int)
    @assert i < j < k "init_loop!: need i < j < k"
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
    b.added_edgeset     = BitSet128()          # test version: edgesets left empty
    b.forbidden_edgeset = BitSet128()          # test version: edgesets left empty
    return b
end

# ── init_loop! (search version) ───────────────────────────────────────────────
# Takes a renumbered triangle index t and precomputed edgesets from the main
# search loop. Use this form inside backtrack!.
function init_loop!(b       ::BdryLoop,
                    t       ::Int16,
                    tmap    ::Vector{NTuple{3,Int}},
                    edgesets::Vector{Tri_Edgesets})
    @inbounds begin
        fill!(b.status, UNUSED)
        i,j,k = tmap[t]
        b.v[i] = BdryVE(Int8(j), Int8(k), t)
        b.v[j] = BdryVE(Int8(k), Int8(i), t)
        b.v[k] = BdryVE(Int8(i), Int8(j), t)
        b.status[i] = ON_BOUNDARY
        b.status[j] = ON_BOUNDARY
        b.status[k] = ON_BOUNDARY
    end
    b.head              = i
    b.n_unused          = b.n - 3
    b.added_edgeset     = edgesets[t].has
    b.forbidden_edgeset = edgesets[t].conf
    return b
end
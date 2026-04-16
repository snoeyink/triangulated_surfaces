# ─────────────────────────────────────────────────────────────────────────────
# BdryLoop.jl
# External dependencies (must be included/loaded before this file):
#   BitSet128.jl   — provides BitSet128
#   triangle_index(i,j,k)::UInt16  requires i < j < k
# ─────────────────────────────────────────────────────────────────────────────

# ── BdryVE ───────────────────────────────────────────────────────────────────
# sizeof == 4, isbits, no padding.
# Vertex i owns directed boundary edge i → next.
# Pending left triangle is (i, next, opp) with renumbered index tri.
const LINK = UInt8(0) # link if nod ear0 or ear1
const EAR0 = UInt8(2) # ear if opp(i) == nxt2(i)
const EAR1 = UInt8(3) # ear if nxt(opp(i)) == i

struct BdryVE
    tri ::UInt16
    next::UInt8
    opp ::UInt8
    type::UInt8
end

# ── Vertex status ─────────────────────────────────────────────────────────────
const UNUSED      = UInt16(0)
const INTERIOR    = UInt16(1024)
const ON_BDRY_MASK = UInt16(0x3FF) # 10 bits for tmax or min triangle index for non-removabble links with this opp

struct VState
    status::UInt16 # status, with min rem link having this opp 
end
zeros(::Type{VState}, n::Integer) = VState.(zeros(UInt16, n))

# ── BdryLoop ──────────────────────────────────────────────────────────────────
struct BdryLoop
    ve   ::Vector{BdryVE}  # half-edge data, length n
    vs   ::Vector{VState}  # UNUSED / ON_BOUNDARY / INTERIOR
    ti :: TriTable
    vunused::BitSet16
end

function BdryLoop(n::Integer)
    BdryLoop(
        Vector{BdryVE}(undef, n),
        zeros(VState, n),
        TriTable(),
        BitSet16(1<<n-1) # max n == 16
end

# ── Accessors ────────────────────────────────────────────────────────────────
@inline nxt(b::BdryLoop, i::Integer)  = @inbounds b.ve[i].next
@inline opp(b::BdryLoop, i::Integer)  = @inbounds b.ve[i].opp
@inline tri(b::BdryLoop, i::Integer)  = @inbounds b.ve[i].tri
@inline nxt2(b::BdryLoop, i::Integer) = nxt(b, nxt(b, i))
@inline nxt3(b::BdryLoop, i::Integer) = nxt(b, nxt2(b, i))

# ── Status predicates ─────────────────────────────────────────────────────────
@inline unused(b::BdryLoop, i::Integer)      = @inbounds b.vs[i].status == UNUSED
@inline on_bdry(b::BdryLoop, i::Integer) = @inbounds b.vs[i].status & ON_BDRY_MASK != 0
@inline interior(b::BdryLoop, i::Integer)    = @inbounds b.vs[i].status == INTERIOR

# ── Boundary vertex predicates ────────────────────────────────────────────────
@inline ear0(b::BdryLoop, i::Integer) = opp(b, i) == nxt2(b, i)
@inline ear1(b::BdryLoop, i::Integer) = nxt(b, opp(b, i)) == i
@inline link(b::BdryLoop, i::Integer) = !ear0(b, i) && !ear1(b, i)

# Removable: safe to add the pending triangle (no non-manifold edge).
# Not removable iff link(i) AND opp(i) is already on the boundary.
@inline removable(b::BdryLoop, i::Integer) = !link(b, i) || !on_boundary(b, opp(b, i))

# first_tri: true for exactly one vertex when the boundary is a single triangle.
@inline function first_tri(b::BdryLoop, i::Integer)
    j = nxt(b, i); k = nxt(b, j)
    opp(b, i) == k && nxt(b, k) == i && i < j && i < k
end

# ── init_loop! (test version) ─────────────────────────────────────────────────
# Uses triangle_index directly; leaves edgesets zeroed.
# Call this form from tests and anywhere renumbering is not needed.
function init_loop!(b::BdryLoop, i::Integer, j::Integer, k::Integer)
    @assert i < j < k "init_loop!: need i < j < k"
    t = triangle_index(i, j, k)
    @inbounds begin
        b.ve[i] = BdryVE(UInt8(j), UInt8(k), t)
        b.ve[j] = BdryVE(UInt8(k), UInt8(i), t)
        b.ve[k] = BdryVE(UInt8(i), UInt8(j), t)
        b.vs[i] = VState(t)
        b.vs[j] = VState(t)
        b.vs[k] = VState(t)
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
                    t       ::UInt16,
                    tmap    ::Vector{TriIJK},
                    edgesets::Vector{Tri_Edgesets})
    @inbounds begin
        fill!(b.vs, VState(UNUSED))
        i,j,k = tmap[t]
        b.ve[i] = BdryVE(j, k, t)
        b.ve[j] = BdryVE(k, i, t)
        b.ve[k] = BdryVE(i, j, t)
        b.vs[i] = VState(t)
        b.vs[j] = VState(t)
        b.vs[k] = VState(t)
    end
    b.head              = i
    b.n_unused          = b.n - 3
    b.added_edgeset     = edgesets[t].has
    b.forbidden_edgeset = edgesets[t].conf
    return b, es
end
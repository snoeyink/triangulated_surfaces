# test/bdry_loop.jl
# Tests for BdryLoop.jl only.
# No BdryStack, no Backtrack.jl dependencies.
# make_six_vertex_loop is defined here and used by test/backtrack.jl
# via include order in runtests.jl.

const BdryVE       = TriangulatedSurfaces.BdryVE
const BdryLoop     = TriangulatedSurfaces.BdryLoop
const UNUSED       = TriangulatedSurfaces.UNUSED
const ON_BOUNDARY  = TriangulatedSurfaces.ON_BOUNDARY
const INTERIOR     = TriangulatedSurfaces.INTERIOR
const MAX_VERTICES = TriangulatedSurfaces.MAX_VERTICES

const triangle_index = TriangulatedSurfaces.triangle_index
const nxt          = TriangulatedSurfaces.nxt
const opp          = TriangulatedSurfaces.opp
const tri          = TriangulatedSurfaces.tri
const nxt2         = TriangulatedSurfaces.nxt2
const nxt3         = TriangulatedSurfaces.nxt3
const unused       = TriangulatedSurfaces.unused
const on_boundary  = TriangulatedSurfaces.on_boundary
const interior     = TriangulatedSurfaces.interior
const ear0         = TriangulatedSurfaces.ear0
const ear1         = TriangulatedSurfaces.ear1
const link         = TriangulatedSurfaces.link
const removable    = TriangulatedSurfaces.removable
const first_tri    = TriangulatedSurfaces.first_tri
const init_loop!   = TriangulatedSurfaces.init_loop!

# ── Helper ────────────────────────────────────────────────────────────────────
# Six-vertex boundary loop 1→3→6→4→7→5→1, interior vertex 2.
# Manually constructed (bypasses init_loop!, so n_unused is not set).
# Used by backtrack.jl tests via include order.
function make_six_vertex_loop(; debug=true)
    b = BdryLoop(16; debug)
    @inbounds begin
        b.v[1] = BdryVE(Int8(3), Int8(5), triangle_index(1,3,5))
        b.v[3] = BdryVE(Int8(6), Int8(2), triangle_index(2,3,6))
        b.v[6] = BdryVE(Int8(4), Int8(7), triangle_index(4,6,7))
        b.v[4] = BdryVE(Int8(7), Int8(6), triangle_index(4,6,7))
        b.v[7] = BdryVE(Int8(5), Int8(6), triangle_index(5,6,7))
        b.v[5] = BdryVE(Int8(1), Int8(3), triangle_index(1,3,5))
        b.v[2] = BdryVE(Int8(6), Int8(5), triangle_index(2,5,6))   # preserved interior state
        for i in (1,3,4,5,6,7); b.status[i] = ON_BOUNDARY; end
        b.status[2] = INTERIOR
    end
    b.head = 1
    return b
end

@testset "BdryLoop" begin

# ─── layout ──────────────────────────────────────────────────────────────────
@testset "BdryVE isbits and size" begin
    @test isbitstype(BdryVE)
    @test sizeof(BdryVE) == 4           # Int8 + Int8 + Int16, no padding
end

# ─── init_loop!(b, 1, 2, 3) ──────────────────────────────────────────────────
@testset "init_loop!(b, 1, 2, 3)" begin
    b = BdryLoop(8)
    init_loop!(b, 1, 2, 3)
    t = triangle_index(1, 2, 3)

    # Half-edge pointers
    @test nxt(b,1)==2  && nxt(b,2)==3  && nxt(b,3)==1
    @test opp(b,1)==3  && opp(b,2)==1  && opp(b,3)==2
    @test tri(b,1)==t  && tri(b,2)==t  && tri(b,3)==t
    @test nxt2(b,1)==3 && nxt2(b,2)==1 && nxt2(b,3)==2
    @test nxt3(b,1)==1 && nxt3(b,2)==2 && nxt3(b,3)==3

    # Status
    @test on_boundary(b,1) && on_boundary(b,2) && on_boundary(b,3)
    @test all(unused(b,i) for i in 4:b.n)

    # New fields
    @test b.n        == 8
    @test b.n_unused == 5              # n - 3 = 8 - 3
    @test iszero(b.added_edgeset)      # test version: edgesets left empty
    @test iszero(b.forbidden_edgeset)
    
    # Single triangle: every vertex is ear0 ∧ ear1 ∧ ¬link ∧ removable
    for i in 1:3
        @test ear0(b,i) && ear1(b,i) && !link(b,i) && removable(b,i)
    end

    # first_tri: only the minimum-indexed vertex qualifies
    @test  first_tri(b,1)
    @test !first_tri(b,2)              # 2 < opp(2)=1 fails
    @test !first_tri(b,3)              # 3 < nxt(3)=1 fails

    @test b.head == 1

    @test 0 == @allocated init_loop!(b, 1, 2, 3)
end

# ─── six-vertex loop predicates ───────────────────────────────────────────────
# Boundary 1→3→6→4→7→5→1, interior vertex 2.
# Verified by hand:
#   ear0: 6 (opp=7=nxt2(6)=5... wait let me re-check)
#   v[6]: next=4, opp=7; nxt2(6)=nxt(4)=7 → opp==nxt2 → ear0(6) ✓
#   v[5]: next=1, opp=3; nxt2(5)=nxt(1)=3 → opp==nxt2 → ear0(5) ✓
#   ear1: 1 (nxt(opp(1))=nxt(5)=1 ✓), 4 (nxt(opp(4))=nxt(6)=4 ✓)
#   link: 3 (opp=2 interior → removable), 7 (opp=6 on boundary → NOT removable)
@testset "six-vertex loop: predicates" begin
    b = make_six_vertex_loop()

    @test !ear0(b,1) && !ear0(b,3) &&  ear0(b,6)
    @test !ear0(b,4) && !ear0(b,7) &&  ear0(b,5)

    @test  ear1(b,1) && !ear1(b,3) && !ear1(b,6)
    @test  ear1(b,4) && !ear1(b,7) && !ear1(b,5)

    @test !link(b,1) &&  link(b,3) && !link(b,6)
    @test !link(b,4) &&  link(b,7) && !link(b,5)

    @test  removable(b,3)    # link, opp=2 is interior
    @test !removable(b,7)    # link, opp=6 is on boundary
    @test  removable(b,1) && removable(b,6) && removable(b,4) && removable(b,5)

    # Not a triangle, not 3-vertex loop
    for i in (1,3,4,5,6,7)
        @test !first_tri(b,i) && nxt3(b,i) != i
    end
end

end # @testset "BdryLoop"
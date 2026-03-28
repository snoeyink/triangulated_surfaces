# test/bdry_loop.jl
const BdryVE      = TriangulatedSurfaces.BdryVE
const BdryLoop    = TriangulatedSurfaces.BdryLoop
const BdryStack   = TriangulatedSurfaces.BdryStack
const UNUSED      = TriangulatedSurfaces.UNUSED
const ON_BOUNDARY = TriangulatedSurfaces.ON_BOUNDARY
const INTERIOR    = TriangulatedSurfaces.INTERIOR
const MAX_VERTICES = TriangulatedSurfaces.MAX_VERTICES

const triangle_index = TriangulatedSurfaces.triangle_index
const t_index    = TriangulatedSurfaces.t_index
const nxt        = TriangulatedSurfaces.nxt
const opp        = TriangulatedSurfaces.opp
const tri        = TriangulatedSurfaces.tri
const nxt2       = TriangulatedSurfaces.nxt2
const nxt3       = TriangulatedSurfaces.nxt3
const unused     = TriangulatedSurfaces.unused
const on_boundary = TriangulatedSurfaces.on_boundary
const interior   = TriangulatedSurfaces.interior
const ear0       = TriangulatedSurfaces.ear0
const ear1       = TriangulatedSurfaces.ear1
const link       = TriangulatedSurfaces.link
const removable  = TriangulatedSurfaces.removable
const first_tri  = TriangulatedSurfaces.first_tri
const init_loop! = TriangulatedSurfaces.init_loop!
const add_ear!   = TriangulatedSurfaces.add_ear!
const add_link!  = TriangulatedSurfaces.add_link!
const popBE!     = TriangulatedSurfaces.popBE!

# ─────────────────────────────────────────────────────────────────────────────
# Helper: six-vertex boundary loop 1→3→6→4→7→5→1, interior vertex 2.
# ─────────────────────────────────────────────────────────────────────────────
function make_six_vertex_loop(; debug=true)
    b = BdryLoop(; debug)
    @inbounds begin
        b.v[1] = BdryVE(Int8(3), Int8(5), t_index(1,3,5))
        b.v[3] = BdryVE(Int8(6), Int8(2), t_index(2,3,6))
        b.v[6] = BdryVE(Int8(4), Int8(7), t_index(4,6,7))
        b.v[4] = BdryVE(Int8(7), Int8(6), t_index(4,6,7))
        b.v[7] = BdryVE(Int8(5), Int8(6), t_index(5,6,7))
        b.v[5] = BdryVE(Int8(1), Int8(3), t_index(1,3,5))
        b.v[2] = BdryVE(Int8(6), Int8(5), t_index(2,5,6))   # preserved interior state
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
    @test sizeof(BdryVE) == 4        # Int8 + Int8 + Int16, no padding
end

# ─── init_loop!(b, 1, 2, 3) ──────────────────────────────────────────────────
@testset "init_loop!(b, 1, 2, 3)" begin
    b = BdryLoop()
    init_loop!(b, 1, 2, 3)
    t = triangle_index(1, 2, 3)

    @test nxt(b,1)==2  && nxt(b,2)==3  && nxt(b,3)==1
    @test opp(b,1)==3  && opp(b,2)==1  && opp(b,3)==2
    @test tri(b,1)==t  && tri(b,2)==t  && tri(b,3)==t
    @test nxt2(b,1)==3 && nxt2(b,2)==1 && nxt2(b,3)==2
    @test nxt3(b,1)==1 && nxt3(b,2)==2 && nxt3(b,3)==3

    @test on_boundary(b,1) && on_boundary(b,2) && on_boundary(b,3)
    @test all(unused(b,i) for i in 4:MAX_VERTICES)

    for i in 1:3
        @test ear0(b,i) && ear1(b,i) && !link(b,i) && removable(b,i)
    end

    @test  first_tri(b,1)
    @test !first_tri(b,2)     # 2 < nxt2(2)=1 fails
    @test !first_tri(b,3)     # 3 < nxt(3)=1 fails

    @test b.head == 1

    @test 0 == @allocated init_loop!(b, 1, 2, 3)
end

# ─── add_ear!(b, s, 3, 4) ────────────────────────────────────────────────────
# j = nxt(3) = 1,  k = 4,  triangle (1,3,4)
# Boundary becomes 1→2→3→4→1.
@testset "add_ear!(b, s, 3, 4)" begin
    b = BdryLoop(); s = BdryStack()
    init_loop!(b, 1, 2, 3)
    add_ear!(b, s, 3, 4)
    t = t_index(1, 3, 4)

    @test nxt(b,3)==4  && opp(b,3)==1  && tri(b,3)==t
    @test nxt(b,4)==1  && opp(b,4)==3  && tri(b,4)==t
    @test on_boundary(b,4) && s.sp==1

    @test  ear0(b,3) && !ear1(b,3)
    @test !ear0(b,4) &&  ear1(b,4)
    @test  removable(b,3) && removable(b,4)

    @test nxt(b,1)==2 && opp(b,1)==3   # v[1] untouched
    @test nxt(b,2)==3 && opp(b,2)==1   # v[2] untouched

    @test nxt(b,1)==2&&nxt(b,2)==3&&nxt(b,3)==4&&nxt(b,4)==1

    @test nxt3(b,1) != 1               # no longer a 3-vertex loop

    b2 = BdryLoop(debug=false); s2 = BdryStack(debug=false)
    init_loop!(b2, 1, 2, 3)
    @test 0 == @allocated add_ear!(b2, s2, 3, 4)
end

# ─── popBE!(b, s, 3) — ear branch ────────────────────────────────────────────
@testset "popBE!(b, s, 3) — ear branch" begin
    b = BdryLoop(); s = BdryStack()
    init_loop!(b, 1, 2, 3); add_ear!(b, s, 3, 4)
    popBE!(b, s, 3)

    @test nxt(b,3)==1  && opp(b,3)==2
    @test tri(b,3) == triangle_index(1, 2, 3)
    @test unused(b,4) && s.sp==0

    b2 = BdryLoop(debug=false); s2 = BdryStack(debug=false)
    init_loop!(b2, 1, 2, 3); add_ear!(b2, s2, 3, 4)
    @test 0 == @allocated popBE!(b2, s2, 3)
end

# ─── two add_ear! calls ───────────────────────────────────────────────────────
# add_ear!(3,4): boundary 1→2→3→4→1
# add_ear!(4,5): j=nxt(4)=1, k=5; triangle(1,4,5); boundary 1→2→3→4→5→1
@testset "two add_ear! calls: 1→2→3→4→5→1" begin
    b = BdryLoop(); s = BdryStack()
    init_loop!(b, 1, 2, 3)
    add_ear!(b, s, 3, 4)
    add_ear!(b, s, 4, 5)
    t = t_index(1, 4, 5)

    @test nxt(b,4)==5  && opp(b,4)==1  && tri(b,4)==t
    @test nxt(b,5)==1  && opp(b,5)==4  && tri(b,5)==t
    @test on_boundary(b,5) && s.sp==2

    @test  ear0(b,4) && !ear1(b,4) && removable(b,4)
    @test !ear0(b,5) &&  ear1(b,5) && removable(b,5)

    @test nxt(b,1)==2&&nxt(b,2)==3&&nxt(b,3)==4&&nxt(b,4)==5&&nxt(b,5)==1
end

# ─── six-vertex loop: ear0, ear1, link, removable ────────────────────────────
# Boundary 1→3→6→4→7→5→1, interior vertex 2.
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

    for i in (1,3,4,5,6,7)
        @test !first_tri(b,i) && nxt3(b,i) != i
    end
end

# ─── add_link!(b, s, 3) ──────────────────────────────────────────────────────
# j = nxt(3) = 6,  k = nxt2(3) = 4,  triangle (3,4,6)
# Post: boundary 1→3→4→7→5→1
@testset "add_link!(b, s, 3)" begin
    b = make_six_vertex_loop(); s = BdryStack()
    add_link!(b, s, 3)
    t = t_index(3, 4, 6)            # = triangle_index(3,4,6) = 16

    @test nxt(b,3)==4  && opp(b,3)==6  && tri(b,3)==t
    @test interior(b,6) && s.sp==1

    @test link(b,3) && removable(b,3)

    @test nxt(b,6)==4 && opp(b,6)==7   # v[6] preserved untouched

    @test nxt(b,1)==3&&nxt(b,3)==4&&nxt(b,4)==7&&nxt(b,7)==5&&nxt(b,5)==1

    b2 = make_six_vertex_loop(debug=false); s2 = BdryStack(debug=false)
    @test 0 == @allocated add_link!(b2, s2, 3)
end

# ─── popBE!(b, s, 3) — link branch ───────────────────────────────────────────
@testset "popBE!(b, s, 3) — link branch" begin
    b = make_six_vertex_loop(); s = BdryStack()
    add_link!(b, s, 3)
    popBE!(b, s, 3)

    @test nxt(b,3)==6  && opp(b,3)==2
    @test tri(b,3) == t_index(2, 3, 6)    # = triangle_index(2,3,6) = 13
    @test on_boundary(b,6) && s.sp==0

    @test link(b,3) && removable(b,3)     # vertex 2 still interior
    @test ear0(b,6) && !ear1(b,6)         # v[6] unchanged throughout

    b2 = make_six_vertex_loop(debug=false); s2 = BdryStack(debug=false)
    add_link!(b2, s2, 3)
    @test 0 == @allocated popBE!(b2, s2, 3)
end

# ─── debug stack index tracking ───────────────────────────────────────────────
@testset "debug stack index tracking" begin
    b = BdryLoop(debug=true); s = BdryStack(debug=true)
    init_loop!(b, 1, 2, 3)

    add_ear!(b, s, 3, 4)
    @test Int(s.dbg_idx[1])==3 && s.sp==1

    add_ear!(b, s, 4, 5)
    @test Int(s.dbg_idx[2])==4 && s.sp==2

    popBE!(b, s, 4)
    @test unused(b,5) && s.sp==1

    popBE!(b, s, 3)
    @test unused(b,4) && s.sp==0

    @test nxt(b,3)==1 && opp(b,3)==2
    @test tri(b,3) == triangle_index(1, 2, 3)
end

end # @testset "BdryLoop"
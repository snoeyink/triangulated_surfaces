# test/bdry_loop.jl
#
# Dependencies (provided externally):
#   triangle_index(a,b,c)::Int16  requires a < b < c
#   t_index(a,b,c)::Int16         sorts internally
#
# All hand-verified values are re-derived in comments.

# ─────────────────────────────────────────────────────────────────────────────
# Helper: six-vertex boundary loop 1→3→6→4→7→5→1, interior vertex 2.
# Pending triangles stored in each BdryVE.tri:
#   v[1]: (1,3,5)  v[3]: (2,3,6)  v[6]: (4,6,7)
#   v[4]: (4,6,7)  v[7]: (5,6,7)  v[5]: (1,3,5)
#   v[2]: (2,5,6)  ← preserved interior state
# ─────────────────────────────────────────────────────────────────────────────
function make_six_vertex_loop(; debug=true)
    b = BdryLoop(16; debug)
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
    b.sp = 0; b.head = 1
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
    b = BdryLoop(8)
    init_loop!(b, 1, 2, 3)
    t = triangle_index(1, 2, 3)     # = 1 + 0 + 1 = 2

    # half-edge pointers
    @test nxt(b,1)==2  && nxt(b,2)==3  && nxt(b,3)==1
    @test opp(b,1)==3  && opp(b,2)==1  && opp(b,3)==2
    @test tri(b,1)==t  && tri(b,2)==t  && tri(b,3)==t
    @test nxt2(b,1)==3 && nxt2(b,2)==1 && nxt2(b,3)==2
    @test nxt3(b,1)==1 && nxt3(b,2)==2 && nxt3(b,3)==3

    # status
    @test on_boundary(b,1) && on_boundary(b,2) && on_boundary(b,3)
    @test all(unused(b,i) for i in 4:8)

    # single triangle: every vertex is ear0 ∧ ear1 ∧ ¬link ∧ removable
    for i in 1:3
        @test ear0(b,i) && ear1(b,i) && !link(b,i) && removable(b,i)
    end

    # first_tri: only the minimum-indexed vertex qualifies
    @test  first_tri(b,1)
    @test !first_tri(b,2)     # 2 < nxt2(2)=1 fails
    @test !first_tri(b,3)     # 3 < nxt(3)=1 fails

    # is_last: same condition for a 3-vertex loop
    @test  is_last(b,1)
    @test !is_last(b,2)
    @test !is_last(b,3)

    @test b.sp == 0 && b.head == 1

    @test 0 == @allocated init_loop!(b, 1, 2, 3)
end

# ─── add_ear!(b, 3, 4) ───────────────────────────────────────────────────────
# j = nxt(3) = 1,  k = 4,  triangle (1,3,4)
# Boundary becomes 1→2→3→4→1.
@testset "add_ear!(b, 3, 4)" begin
    b = BdryLoop(8); init_loop!(b, 1, 2, 3)
    add_ear!(b, 3, 4)
    t = t_index(1, 3, 4)             # = triangle_index(1,3,4) = 1+1+1 = 3

    @test nxt(b,3)==4  && opp(b,3)==1  && tri(b,3)==t
    @test nxt(b,4)==1  && opp(b,4)==3  && tri(b,4)==t
    @test on_boundary(b,4) && b.sp==1

    # spec post-conditions
    @test  ear0(b,3) && !ear1(b,3)   # opp(3)=1=nxt2(3); nxt(1)=2≠3
    @test !ear0(b,4) &&  ear1(b,4)   # opp(4)=3; nxt(3)=4=4 ✓
    @test  removable(b,3) && removable(b,4)

    # vertices untouched by add_ear!(3,4)
    @test nxt(b,1)==2 && opp(b,1)==3
    @test nxt(b,2)==3 && opp(b,2)==1

    # full boundary traversal: 1→2→3→4→1
    @test nxt(b,1)==2&&nxt(b,2)==3&&nxt(b,3)==4&&nxt(b,4)==1

    # no longer a single triangle
    @test !is_last(b,1) && !first_tri(b,1)

    b2 = BdryLoop(8, debug=false); init_loop!(b2, 1, 2, 3)
    @test 0 == @allocated add_ear!(b2, 3, 4)
end

# ─── popBE!(b, 3) — ear branch ───────────────────────────────────────────────
# ear0(3) = true after add_ear!, so popBE! marks nxt(3)=4 as unused.
@testset "popBE!(b, 3) — ear branch" begin
    b = BdryLoop(8); init_loop!(b, 1, 2, 3); add_ear!(b, 3, 4)
    popBE!(b, 3)

    @test nxt(b,3)==1  && opp(b,3)==2            # v[3] restored to post-init state
    @test tri(b,3) == triangle_index(1, 2, 3)
    @test unused(b,4) && b.sp==0

    b2 = BdryLoop(8, debug=false); init_loop!(b2, 1, 2, 3); add_ear!(b2, 3, 4)
    @test 0 == @allocated popBE!(b2, 3)
end

# ─── two add_ear! calls ───────────────────────────────────────────────────────
# add_ear!(3,4): boundary 1→2→3→4→1
# add_ear!(4,5): j=nxt(4)=1, k=5; triangle(1,4,5); boundary 1→2→3→4→5→1
@testset "two add_ear! calls: 1→2→3→4→5→1" begin
    b = BdryLoop(8); init_loop!(b, 1, 2, 3)
    add_ear!(b, 3, 4)
    add_ear!(b, 4, 5)
    t = t_index(1, 4, 5)

    @test nxt(b,4)==5  && opp(b,4)==1  && tri(b,4)==t
    @test nxt(b,5)==1  && opp(b,5)==4  && tri(b,5)==t
    @test on_boundary(b,5) && b.sp==2

    @test  ear0(b,4) && !ear1(b,4) && removable(b,4)
    @test !ear0(b,5) &&  ear1(b,5) && removable(b,5)

    # full boundary traversal
    @test nxt(b,1)==2&&nxt(b,2)==3&&nxt(b,3)==4&&nxt(b,4)==5&&nxt(b,5)==1
end

# ─── six-vertex loop: ear0, ear1, link, removable ────────────────────────────
# Boundary 1→3→6→4→7→5→1, interior vertex 2.
# Verified by hand:
#   ear0: 6 (opp=7=nxt2),  5 (opp=3=nxt2)
#   ear1: 1 (nxt(opp=5)=1), 4 (nxt(opp=6)=4)
#   link: 3 (opp=2 interior → removable)
#         7 (opp=6 on boundary → NOT removable)
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

    # six-vertex loop: no first_tri or is_last
    for i in (1,3,4,5,6,7)
        @test !first_tri(b,i) && !is_last(b,i)
    end
end

# ─── add_link!(b, 3) ─────────────────────────────────────────────────────────
# j = nxt(3) = 6,  k = nxt2(3) = 4,  triangle (3,4,6)
# SPEC CORRECTION: v[i] = BdryVE(k, j, t), NOT BdryVE(j, k, t).
#   Wrong form BdryVE(j=6, k=4, t): nxt=6(interior!), opp=4;
#     nxt2(3)=nxt(6)=4=opp → ear0(3)=true, killing link post-condition
#     and sending popBE! into the wrong branch.
#   Correct form BdryVE(k=4, j=6, t): nxt=4(boundary), opp=6(interior);
#     ear0(3): opp=6 ≠ nxt2(3)=7 → false ✓
#     ear1(3): nxt(6)=4 ≠ 3 → false ✓  → link(3)=true ✓
# Post: boundary 1→3→4→7→5→1
@testset "add_link!(b, 3)" begin
    b = make_six_vertex_loop()
    add_link!(b, 3)
    t = t_index(3, 4, 6)            # = triangle_index(3,4,6) = 3+3+10 = 16

    @test nxt(b,3)==4  && opp(b,3)==6  && tri(b,3)==t
    @test interior(b,6) && b.sp==1

    # post-condition: link(3) still holds
    @test link(b,3)
    @test removable(b,3)             # opp=6 is now interior

    # v[6] preserved untouched (critical for undo)
    @test nxt(b,6)==4 && opp(b,6)==7

    # boundary traversal: 1→3→4→7→5→1
    @test nxt(b,1)==3&&nxt(b,3)==4&&nxt(b,4)==7&&nxt(b,7)==5&&nxt(b,5)==1

    b2 = make_six_vertex_loop(debug=false)
    @test 0 == @allocated add_link!(b2, 3)
end

# ─── popBE!(b, 3) — link branch ──────────────────────────────────────────────
# link(3)=true, opp(3)=6 interior → link branch:
#   mark 6 on_boundary; restore v[3] from stack.
@testset "popBE!(b, 3) — link branch" begin
    b = make_six_vertex_loop()
    add_link!(b, 3)
    popBE!(b, 3)

    @test nxt(b,3)==6  && opp(b,3)==2     # v[3] restored
    @test tri(b,3) == t_index(2, 3, 6)    # = triangle_index(2,3,6) = 13
    @test on_boundary(b,6) && b.sp==0

    # state matches pre-add_link! exactly
    @test link(b,3) && removable(b,3)     # vertex 2 still interior
    @test ear0(b,6) && !ear1(b,6)         # v[6] unchanged throughout

    b2 = make_six_vertex_loop(debug=false); add_link!(b2, 3)
    @test 0 == @allocated popBE!(b2, 3)
end

# ─── debug index tracking ─────────────────────────────────────────────────────
@testset "debug stack index tracking" begin
    b = BdryLoop(8, debug=true)
    init_loop!(b, 1, 2, 3)          # sp=0

    add_ear!(b, 3, 4)                # push for vertex 3; sp→1
    @test Int(b.dbg_idx[1])==3 && b.sp==1

    add_ear!(b, 4, 5)                # push for vertex 4; sp→2; j=1, k=5
    @test Int(b.dbg_idx[2])==4 && b.sp==2

    popBE!(b, 4)                     # ear0(4): nxt(4)=5 → 5 unused
    @test unused(b,5) && b.sp==1

    popBE!(b, 3)                     # ear0(3): nxt(3)=4 → 4 unused
    @test unused(b,4) && b.sp==0

    # fully restored to post-init_loop! state
    @test nxt(b,3)==1 && opp(b,3)==2
    @test tri(b,3) == triangle_index(1, 2, 3)
end

# ─── is_last ──────────────────────────────────────────────────────────────────
@testset "is_last" begin
    b = BdryLoop(5); init_loop!(b, 1, 2, 3)
    @test  is_last(b,1) && !is_last(b,2) && !is_last(b,3)

    b2 = BdryLoop(6); init_loop!(b2, 2, 3, 5)   # non-unit starting indices
    @test  is_last(b2,2) && !is_last(b2,3) && !is_last(b2,5)

    add_ear!(b, 3, 4)   # boundary 1→2→3→4→1: no longer 3-vertex
    @test !is_last(b,1) && !is_last(b,3)
end

end # @testset "BdryLoop"
# test/backtrack.jl
# Tests for Backtrack.jl.
# Requires bdry_loop.jl to be included first in runtests.jl:
#   provides BdryVE, BdryLoop, UNUSED, ON_BOUNDARY, INTERIOR, MAX_VERTICES,
#   triangle_index, t_index, nxt/opp/tri/nxt2/nxt3, status predicates,
#   ear0/ear1/link/removable/first_tri, init_loop!, and make_six_vertex_loop.

const BitSet128        = TriangulatedSurfaces.BitSet128
const BdryStack        = TriangulatedSurfaces.BdryStack
const SearchFrame      = TriangulatedSurfaces.SearchFrame
const add_ear!         = TriangulatedSurfaces.add_ear!
const add_link!        = TriangulatedSurfaces.add_link!
const popBE!           = TriangulatedSurfaces.popBE!
const ear_ok           = TriangulatedSurfaces.ear_ok
const link_ok          = TriangulatedSurfaces.link_ok
const is_min_removable = TriangulatedSurfaces.is_min_removable
const edge_bit         = TriangulatedSurfaces.edge_bit

# ── Test fixtures ─────────────────────────────────────────────────────────────
# Tri-table using original t_index ordering (no renumbering).
# All 6 permutations of each valid triple are filled; degenerate entries stay 0.
function _make_test_tri_table()
    T = zeros(Int16, MAX_VERTICES, MAX_VERTICES, MAX_VERTICES)
    for a in 1:MAX_VERTICES, b in 1:MAX_VERTICES, c in 1:MAX_VERTICES
        if a != b && b != c && a != c
            @inbounds T[a,b,c] = Int16(t_index(a, b, c))
        end
    end
    return T
end

const TEST_TRI   = _make_test_tri_table()
# All-zero conflict edgesets → no geometry conflicts; tests topology only.
const TEST_ESETS = fill(BitSet128(), Int(triangle_index(14, 15, 16)))

@testset "Backtrack" begin

# ─── layout ──────────────────────────────────────────────────────────────────
@testset "SearchFrame isbits and size" begin
    @test isbitstype(SearchFrame)
    @test sizeof(SearchFrame) == 48    # 2×16 + 4 + 5×1 = 41, padded to 48
end

# ─── add_ear!(b, s, 3, 4) ────────────────────────────────────────────────────
# j = nxt(3) = 1, k = 4; triangle (1,3,4); boundary becomes 1→2→3→4→1.
@testset "add_ear!(b, s, 3, 4)" begin
    b = BdryLoop(8); s = BdryStack()
    init_loop!(b, 1, 2, 3)
    add_ear!(b, s, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    t = t_index(1, 3, 4)

    @test nxt(b,3)==4  && opp(b,3)==1  && tri(b,3)==t
    @test nxt(b,4)==1  && opp(b,4)==3  && tri(b,4)==t
    @test on_boundary(b,4) && s.sp==1
    @test b.n_unused == 4              # decremented from 5

    # Stack frame identity
    @test Int(s.frames[1].vertex_i) == 3

    # Spec post-conditions
    @test  ear0(b,3) && !ear1(b,3)
    @test !ear0(b,4) &&  ear1(b,4)
    @test  removable(b,3) && removable(b,4)

    # Untouched vertices
    @test nxt(b,1)==2 && opp(b,1)==3
    @test nxt(b,2)==3 && opp(b,2)==1

    # Full boundary traversal: 1→2→3→4→1
    @test nxt(b,1)==2 && nxt(b,2)==3 && nxt(b,3)==4 && nxt(b,4)==1
    @test nxt3(b,1) != 1               # no longer a 3-vertex loop

    b2 = BdryLoop(8; debug=false); s2 = BdryStack(debug=false)
    init_loop!(b2, 1, 2, 3)
    @test 0 == @allocated add_ear!(b2, s2, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
end

# ─── popBE!(b, s, 3) — ear branch ────────────────────────────────────────────
@testset "popBE!(b, s, 3) — ear branch" begin
    b = BdryLoop(8); s = BdryStack()
    init_loop!(b, 1, 2, 3)
    add_ear!(b, s, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    popBE!(b, s, 3)

    @test nxt(b,3)==1  && opp(b,3)==2
    @test tri(b,3) == triangle_index(1, 2, 3)
    @test unused(b,4) && s.sp==0
    @test b.n_unused == 5              # restored

    b2 = BdryLoop(8; debug=false); s2 = BdryStack(debug=false)
    init_loop!(b2, 1, 2, 3); add_ear!(b2, s2, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    @test 0 == @allocated popBE!(b2, s2, 3)
end

# ─── two add_ear! calls: 1→2→3→4→5→1 ────────────────────────────────────────
# add_ear!(3,4): boundary 1→2→3→4→1
# add_ear!(4,5): j=nxt(4)=1, k=5; triangle (1,4,5); boundary 1→2→3→4→5→1
@testset "two add_ear!: 1→2→3→4→5→1" begin
    b = BdryLoop(8); s = BdryStack()
    init_loop!(b, 1, 2, 3)
    add_ear!(b, s, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    add_ear!(b, s, TEST_TRI, TEST_ESETS, 4, 5, Int8(0), Int8(0))
    t = t_index(1, 4, 5)

    @test nxt(b,4)==5  && opp(b,4)==1  && tri(b,4)==t
    @test nxt(b,5)==1  && opp(b,5)==4  && tri(b,5)==t
    @test on_boundary(b,5) && s.sp==2
    @test b.n_unused == 3

    @test  ear0(b,4) && !ear1(b,4) && removable(b,4)
    @test !ear0(b,5) &&  ear1(b,5) && removable(b,5)

    @test nxt(b,1)==2 && nxt(b,2)==3 && nxt(b,3)==4 && nxt(b,4)==5 && nxt(b,5)==1
end

# ─── add_link!(b, s, 3) ──────────────────────────────────────────────────────
# j = nxt(3) = 6, k = nxt2(3) = 4; triangle (3,4,6)
# Correct form: BdryVE(k,j,t) — see spec correction notes in source.
# Post: boundary 1→3→4→7→5→1
@testset "add_link!(b, s, 3)" begin
    b = make_six_vertex_loop(); s = BdryStack()
    add_link!(b, s, TEST_TRI, TEST_ESETS, 3, Int8(0), Int8(0))
    t = t_index(3, 4, 6)              # = triangle_index(3,4,6) = 16

    @test nxt(b,3)==4  && opp(b,3)==6  && tri(b,3)==t
    @test interior(b,6) && s.sp==1
    @test Int(s.frames[1].vertex_i) == 3

    @test link(b,3) && removable(b,3)

    # v[6] preserved untouched (critical for undo)
    @test nxt(b,6)==4 && opp(b,6)==7

    # Boundary traversal: 1→3→4→7→5→1
    @test nxt(b,1)==3 && nxt(b,3)==4 && nxt(b,4)==7 && nxt(b,7)==5 && nxt(b,5)==1

    b2 = make_six_vertex_loop(debug=false); s2 = BdryStack(debug=false)
    @test 0 == @allocated add_link!(b2, s2, TEST_TRI, TEST_ESETS, 3, Int8(0), Int8(0))
end

# ─── popBE!(b, s, 3) — link branch ───────────────────────────────────────────
@testset "popBE!(b, s, 3) — link branch" begin
    b = make_six_vertex_loop(); s = BdryStack()
    add_link!(b, s, TEST_TRI, TEST_ESETS, 3, Int8(0), Int8(0))
    popBE!(b, s, 3)

    @test nxt(b,3)==6  && opp(b,3)==2
    @test tri(b,3) == t_index(2, 3, 6)    # = triangle_index(2,3,6) = 13
    @test on_boundary(b,6) && s.sp==0

    @test link(b,3) && removable(b,3)     # vertex 2 still interior
    @test ear0(b,6) && !ear1(b,6)         # v[6] unchanged throughout

    b2 = make_six_vertex_loop(debug=false); s2 = BdryStack(debug=false)
    add_link!(b2, s2, TEST_TRI, TEST_ESETS, 3, Int8(0), Int8(0))
    @test 0 == @allocated popBE!(b2, s2, 3)
end

# ─── head update when add_link! buries head ──────────────────────────────────
# Boundary 3→6→4→7→5→3, head=3.
# add_link!(5): j = nxt(5) = 3 = head → head must update to k = nxt2(5) = 6.
@testset "head update when add_link! buries head" begin
    b = BdryLoop(16); s = BdryStack()
    @inbounds begin
        b.v[3] = BdryVE(Int8(6), Int8(5), t_index(3,5,6))
        b.v[6] = BdryVE(Int8(4), Int8(7), t_index(4,6,7))
        b.v[4] = BdryVE(Int8(7), Int8(6), t_index(4,6,7))
        b.v[7] = BdryVE(Int8(5), Int8(6), t_index(5,6,7))
        b.v[5] = BdryVE(Int8(3), Int8(4), t_index(3,4,5))
        for i in (3,4,5,6,7); b.status[i] = ON_BOUNDARY; end
    end
    b.head = 3; b.n_unused = 11

    add_link!(b, s, TEST_TRI, TEST_ESETS, 5, Int8(0), Int8(0))

    @test b.head == 6                  # updated from 3 to k = nxt2(5) = 6
    @test interior(b, 3)               # j = 3 buried

    popBE!(b, s, 5)
    @test b.head == 3                  # saved_head restored
    @test on_boundary(b, 3)
end

# ─── ear_ok / link_ok ────────────────────────────────────────────────────────
@testset "ear_ok / link_ok with empty edgesets" begin
    b = BdryLoop(8); init_loop!(b, 1, 2, 3)

    # Valid triangles pass when no conflicts
    @test  ear_ok(b, TEST_TRI, TEST_ESETS, 3, 4)  # j=1, k=4; tri(1,3,4) ≠ 0
    @test  ear_ok(b, TEST_TRI, TEST_ESETS, 1, 4)  # j=2, k=4; tri(1,2,4) ≠ 0

    # Degenerate: k == i → T[1, nxt(1)=2, 1] = 0 → false
    @test !ear_ok(b, TEST_TRI, TEST_ESETS, 1, 1)

    # forbidden_edgeset blocks edge (1,4)
    # ear_ok(b2, ..., 3, 4): new edges are (j=1,k=4) and (i=3,k=4); (1,4) forbidden → false
    # ear_ok(b2, ..., 1, 4): new edges are (j=2,k=4) and (i=1,k=4); (1,4) forbidden → false
    b2 = BdryLoop(8); init_loop!(b2, 1, 2, 3)
    b2.forbidden_edgeset = edge_bit(1, 4)
    @test !ear_ok(b2, TEST_TRI, TEST_ESETS, 3, 4)
    @test !ear_ok(b2, TEST_TRI, TEST_ESETS, 1, 4)

    # link_ok on six-vertex loop: link(3)=true, new edge (i=3, k=nxt2(3)=4)
    b3 = make_six_vertex_loop()
    @test  link_ok(b3, TEST_TRI, TEST_ESETS, 3)   # edge(3,4) not forbidden

    b3.forbidden_edgeset = edge_bit(3, 4)
    @test !link_ok(b3, TEST_TRI, TEST_ESETS, 3)   # edge(3,4) now forbidden
end

# ─── is_min_removable ────────────────────────────────────────────────────────
@testset "is_min_removable" begin
    b = BdryLoop(8); s = BdryStack()
    init_loop!(b, 1, 2, 3)

    # Single triangle: tri(1)=tri(2)=tri(3)=triangle_index(1,2,3)=2
    # All removable and equal → min satisfied for any
    @test is_min_removable(b, 1)
    @test is_min_removable(b, 2)

    # After add_ear!(3,4): tri(3) = t_index(1,3,4) = 3
    # Vertex 1 and 2 are still removable with tri = 2 < 3
    # → is_min_removable(b, 3) = false
    add_ear!(b, s, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    @test !is_min_removable(b, 3)     # tri(1)=2 < tri(3)=3
    @test  is_min_removable(b, 1)     # tri(1)=2 is minimum
    @test  is_min_removable(b, 2)     # tri(2)=2 is tied for minimum
end

# ─── SearchFrame vertex_i tracking ───────────────────────────────────────────
@testset "SearchFrame vertex_i tracking" begin
    b = BdryLoop(8; debug=true); s = BdryStack(debug=true)
    init_loop!(b, 1, 2, 3)

    add_ear!(b, s, TEST_TRI, TEST_ESETS, 3, 4, Int8(0), Int8(0))
    @test Int(s.frames[1].vertex_i) == 3 && s.sp == 1

    add_ear!(b, s, TEST_TRI, TEST_ESETS, 4, 5, Int8(0), Int8(0))
    @test Int(s.frames[2].vertex_i) == 4 && s.sp == 2

    popBE!(b, s, 4)                    # ear0(4): nxt(4)=5 → 5 unused
    @test unused(b,5) && s.sp == 1

    popBE!(b, s, 3)                    # ear0(3): nxt(3)=4 → 4 unused
    @test unused(b,4) && s.sp == 0

    # Fully restored to post-init_loop! state
    @test nxt(b,3)==1 && opp(b,3)==2
    @test tri(b,3) == triangle_index(1, 2, 3)
    @test b.n_unused == 5
end

end # @testset "Backtrack"
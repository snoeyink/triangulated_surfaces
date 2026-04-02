# test/backtrack.jl
# Tests for Backtrack.jl.
# Requires bdry_loop.jl to be included first in runtests.jl:
#   provides BdryVE, BdryLoop, UNUSED, ON_BOUNDARY, INTERIOR, MAX_VERTICES,
#   triangle_index, triangle_index, nxt/opp/tri/nxt2/nxt3, status predicates,
#   ear0/ear1/link/removable/first_tri, init_loop!, and make_six_vertex_loop.

const BitSet128        = TriangulatedSurfaces.BitSet128
const BdryStack        = TriangulatedSurfaces.BdryStack
const triangle_index   = TriangulatedSurfaces.triangle_index
const SearchFrame      = TriangulatedSurfaces.SearchFrame
const add_ear!         = TriangulatedSurfaces.add_ear!
const add_link!        = TriangulatedSurfaces.add_link!
const popBE!           = TriangulatedSurfaces.popBE!
const ear_ok           = TriangulatedSurfaces.ear_ok
const link_ok          = TriangulatedSurfaces.link_ok
const is_min_removable = TriangulatedSurfaces.is_min_removable
const precompute_conflicts = TriangulatedSurfaces.precompute_conflicts
const build_tri_table = TriangulatedSurfaces.build_tri_table
const Tri_Edgesets     = TriangulatedSurfaces.Tri_Edgesets

# ── Test fixtures ─────────────────────────────────────────────────────────────
# Tri-table using original triangle_index ordering (no renumbering).
# All 6 permutations of each valid triple are filled; degenerate entries stay 0.



@testset "Backtrack" begin

# ─── layout ──────────────────────────────────────────────────────────────────
@testset "SearchFrame isbits and size" begin
    @test isbitstype(SearchFrame)
    @test sizeof(SearchFrame) == 48    # 2×16 + 4 + 5×1 = 41, padded to 48
end

# ─── add_ear!(b, s, 3, 4) ────────────────────────────────────────────────────
# j = nxt(3) = 1, k = 4; triangle (1,3,4); boundary becomes 1→2→3→4→1.
@testset "add_ear!(b, s, 3, 4)" begin
    points = tetrahedron_with_origin(scale=4)
    push!(points, Point3D(2, 1, 0))
    tm,es = precompute_conflicts(points)
    tmax, tmap, esets, tri_table = build_tri_table(Int8(6), Int16(length(es)), tm, es)  # test tri_table and edgesets for add_ear!/add_link!
    @test tmax == tri_table[3, 4, 6]

    b = BdryLoop(6); s = BdryStack()
    init_loop!(b, Int16(tmax), tmap, esets)
    @test nxt(b,3)==4 && opp(b,3)==6 && tri(b,3)==tmax
    @test nxt(b,4)==6 && opp(b,4)==3 && tri(b,4)==tmax
    @test nxt(b,6)==3 && opp(b,6)==4 && tri(b,6)==tmax
    @test nxt3(b,3) == 3  
                 # 3-vertex loop: nxt2(3)=6, nxt3(3)=3
    @test 0 == @allocated add_ear!(b, s, tri_table, esets, 3, 1, Int8(0), Int8(0))
    
    t = tri_table[3, 1, 4]
    @test nxt(b,3)==1  && opp(b,3)==4  && tri(b,3)==t
    @test nxt(b,1)==4  && opp(b,1)==3  && tri(b,1)==t
    @test on_boundary(b,1) && s.sp==1
    @test b.n_unused == 2              # decremented from 3

    # Stack frame identity
    @test Int(s.frames[1].vertex_i) == 3

    # Spec post-conditions
    @test  ear0(b,3) && !ear1(b,3)
    @test !ear0(b,1) &&  ear1(b,1)
    @test  removable(b,3) && removable(b,4)

    # Untouched vertices
    @test nxt(b,4)==6 && opp(b,4)==3
    @test nxt(b,6)==3 && opp(b,6)==4

    # Full boundary traversal: 1→2→3→4→1
    @test nxt(b,3)==1 && nxt(b,1)==4 && nxt(b,4)==6 && nxt(b,6)==3
    @test nxt3(b,3) != 3               # no longer a 3-vertex loop

    @test 0 == @allocated popBE!(b, s, 3)

    @test nxt(b,3)==4 && opp(b,3)==6 && tri(b,3)==tmax
    @test nxt(b,4)==6 && opp(b,4)==3 && tri(b,4)==tmax
    @test nxt(b,6)==3 && opp(b,6)==4 && tri(b,6)==tmax
    @test nxt3(b,3) == 3  

    @test unused(b,1) && s.sp==0
    @test b.n_unused == 3              # restored
end




end # @testset "Backtrack"
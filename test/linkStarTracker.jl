using Test

include("../src/LinkStarTracker.jl")

@testset "LinkStarTracker Implementation" begin
    using .LinkStarTracker

    # Mock edge_index function for N=15
    @inline function mock_edge_index(a::UInt8, b::UInt8)
        return (a < b) ? 2*a + (b-1)*(b-2) - one(UInt8) : 2*b + (a-1)*(a-2)
    end
    @inline mock_edge_index(a::Integer, b::Integer) = mock_edge_index(UInt8(a), UInt8(b))
    
    N = UInt8(15)
    max_idx = mock_edge_index(N, N - one(UInt8))
    tracker = LinkTracker(N, max_idx)

    @testset "Open Chains and Valid Closure" begin
        # Add T1 = (1, 2, 3)
        success, ci, cj, ck = add_triangle!(tracker, 1, 2, 3, mock_edge_index)
        @test success == true
        @test ci == false # Star 1 is open
        
        # Add T2 = (1, 3, 4) -> Connects to T1 at vertex 1
        success, ci, cj, ck = add_triangle!(tracker, 1, 3, 4, mock_edge_index)
        @test success == true
        @test tracker.degree[1] == 2
        
        # Add T3 = (1, 4, 2) -> Closes the loop at vertex 1!
        success, ci, cj, ck = add_triangle!(tracker, 1, 4, 2, mock_edge_index)
        @test success == true
        @test ci == true # Vertex 1 forms a valid loop of size 3!
        @test tracker.degree[1] == 3
        
        # Undo T3
        undo!(tracker, 1, 4, 2, mock_edge_index)
        @test tracker.degree[1] == 2
        @test tracker.star_next[mock_edge_index(1, 4)] == 0 # Lazy deletion works
    end

    @testset "Premature Cycle Pruning" begin
        tracker2 = LinkTracker(N, max_idx)
        
        # T1 = (5, 6, 7)
        add_triangle!(tracker2, 5, 6, 7, mock_edge_index)
        # T2 = (5, 7, 8)
        add_triangle!(tracker2, 5, 7, 8, mock_edge_index)
        
        # Force a premature cycle by closing the loop at 5, but simulating 
        # that 5 was supposed to have more triangles (degree will be 3, 
        # but what if it connects leaving another piece open?)
        # Wait, if we add (5, 8, 6), degree becomes 3, loop is 3 -> Valid.
        # Let's create two disconnected components at vertex 5!
        
        # Comp 1: (5, 6, 7) and (5, 7, 6) -> Closes immediately
        success, ci, cj, ck = add_triangle!(tracker2, 5, 7, 6, mock_edge_index)
        @test success == false # Premature cycle! Size 2 loop, but degree is 3 (because of T2)
        
        undo!(tracker2, 5, 7, 6, mock_edge_index)
        @test tracker2.degree[5] == 2
    end
end
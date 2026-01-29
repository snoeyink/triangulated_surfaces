# TriangulatedSurfaces

Small Julia project to define 3D points and enumerate triangulations of a convex polygonal surface. Assumes points are ordered on the boundary of a convex polygon (all share the same plane).

## Usage

1. Start Julia in this folder: `julia --project=.`
2. Run the example:
   ```julia
   include("example.jl")
   ```
3. Adjust `points` in `example.jl` to explore other convex polygons.

The code returns all triangulations (Catalan many) using vertex 1 as a fan root and combining sub-triangulations recursively.

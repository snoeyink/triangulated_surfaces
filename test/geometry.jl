@testset "Point geometry" begin
    points = tetrahedron_with_origin(scale=4)

    @test points == [
        Point3D(4, 4, 4),
        Point3D(4, -4, -4),
        Point3D(-4, 4, -4),
        Point3D(-4, -4, 4),
        Point3D(0, 0, 0),
    ]

    @test !point_in_tetra(Point3D(-2, -2, -2), points[1:4]...)
    @test point_in_tetra(Point3D(-1, -1, -1), points[1:4]...)
    @test point_in_tetra(Point3D(2, 1, 0), points[1:4]...)

    push!(points, Point3D(-1, -1, -1))
    @test has_coplanar_quad(points)

    pop!(points)
    push!(points, Point3D(1, 0, 0))
    @test has_coplanar_quad(points)

    pop!(points)
    push!(points, Point3D(2, 1, 0))
    @test !has_coplanar_quad(points)
end

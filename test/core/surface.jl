using TensorNetworkCodes
using Test

@testset "surface_code" begin
    surface = surface_code(3)
    @test verify_code(surface)
    @test find_distance_logicals(surface)[1] == 4
end

@testset "rotated surface_code" begin
    rot_surface = rotated_surface_code(3)
    @test verify_code(rot_surface)
    @test find_distance_logicals(rot_surface)[1] == 3
end

using TensorNetworkCodes
using Test

@testset "surface_code" begin
    surface = surface_code(3)
    @test verify(surface)
    @test find_distance_logicals(surface)[1] == 4
end

@testset "almost rotated surface_code" begin
    almost_rot = almost_rotated_surface_code(3,five_qubit_code(),[2,2])
    @test verify(almost_rot)
    @test find_distance_logicals(almost_rot)[1] == 3
end

@testset "rotated surface_code" begin
    rot_surface = rotated_surface_code(3)
    @test verify(rot_surface)
    @test find_distance_logicals(rot_surface)[1] == 3
end

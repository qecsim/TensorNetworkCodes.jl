using TensorNetworkCodes
using Test

@testset "surface_code" begin
    surface = surface_code(3)
    @test verify_code(surface)
    @test find_distance_logicals(surface)[1] == 4
end

using Test
using TensorNetworkCodes

@testset "hello" begin
    @test hello("fella") == "Hello fella"
end

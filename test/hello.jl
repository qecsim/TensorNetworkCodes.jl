using Test
using TensorNetworkCodes
using Qecsim.BasicModels:FiveQubitCode

@testset "hello" begin
    @test hello("fella") == "Hello fella"
    @test hello(FiveQubitCode()) == "Hello 5-qubit code"
end

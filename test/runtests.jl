using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    @safetestset "hello.jl" begin include("hello.jl") end
end

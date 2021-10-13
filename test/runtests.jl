using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    # SimpleCode tests
    @safetestset "simple/functions.jl" begin include("simple/functions.jl") end
    @safetestset "simple/functions_advanced.jl" begin
        include("simple/functions_advanced.jl")
    end
    # Core tests
    @safetestset "core/types.jl" begin include("core/types.jl") end
    @safetestset "core/contract.jl" begin include("core/contract.jl") end
    @safetestset "core/surface.jl" begin include("core/surface.jl") end
end

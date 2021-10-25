using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    @safetestset "pauli_functions.jl" begin include("pauli_functions.jl") end
    @safetestset "code_functions.jl" begin include("code_functions.jl") end
    # Core tests
    @safetestset "core/types.jl" begin include("core/types.jl") end
    @safetestset "core/contract.jl" begin include("core/contract.jl") end
    @safetestset "core/surface.jl" begin include("core/surface.jl") end
end

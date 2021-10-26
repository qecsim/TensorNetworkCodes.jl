using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    @safetestset "types.jl" begin include("types.jl") end
    @safetestset "pauli_functions.jl" begin include("pauli_functions.jl") end
    @safetestset "code_functions.jl" begin include("code_functions.jl") end
    @safetestset "code_graph_functions.jl" begin include("code_graph_functions.jl") end
    # Simple tests
    # Core tests
    @safetestset "core/contract.jl" begin include("core/contract.jl") end
    @safetestset "core/surface.jl" begin include("core/surface.jl") end
end

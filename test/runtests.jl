using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    @safetestset "types.jl" begin include("types.jl") end
    @safetestset "pauli_functions.jl" begin include("pauli_functions.jl") end
    @safetestset "code_functions.jl" begin include("code_functions.jl") end
    @safetestset "code_graph_functions.jl" begin include("code_graph_functions.jl") end
    @safetestset "contraction_functions.jl" begin include("contraction_functions.jl") end
    @safetestset "qecsim_adaptors.jl" begin include("qecsim_adaptors.jl") end
    # Simple tests
    # Core tests
    @safetestset "core/surface.jl" begin include("core/surface.jl") end
end

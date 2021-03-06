using Test
using SafeTestsets

@testset verbose=true "TensorNetworkCodes.jl" begin
    @safetestset "types.jl" begin include("types.jl") end
    @safetestset "pauli_functions.jl" begin include("pauli_functions.jl") end
    @safetestset "code_functions.jl" begin include("code_functions.jl") end
    @safetestset "code_graph_functions.jl" begin include("code_graph_functions.jl") end
    @safetestset "contraction_functions.jl" begin include("contraction_functions.jl") end
    @safetestset "plotting_functions.jl" begin include("plotting_functions.jl") end
    @safetestset "examples/simple.jl" begin include("examples/simple.jl") end
    @safetestset "examples/surface.jl" begin include("examples/surface.jl") end
    @safetestset "itensors_functions.jl" begin include("itensors_functions.jl") end
    @safetestset "TNDecode.jl" begin include("TNDecode.jl") end
    @safetestset "TNDistance.jl" begin include("TNDistance.jl") end
    @safetestset "QecsimAdaptors.jl" begin include("QecsimAdaptors.jl") end
end

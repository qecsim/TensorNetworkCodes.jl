using TensorNetworkCodes
using Test

@testset "plot_code" begin
    code = TensorNetworkCode(five_qubit_code())
    plot_code(code)  # no error thrown
end

@testset "plot_operator" begin
    code = TensorNetworkCode(steane_code())
    plot_operator(code, code.stabilizers[1])  # no error thrown
end

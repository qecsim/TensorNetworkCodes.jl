using TensorNetworkCodes
using Test

@testset "TensorNetworkCode" begin
    # test TN_code_types.jl
    code = five_qubit_code()
    TN_code1 = TensorNetworkCode(code)
    @test verify_code(TN_code1)
end

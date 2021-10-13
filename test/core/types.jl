using TensorNetworkCodes
using Test

@testset "TNCode" begin
    # test TN_code_types.jl
    code = five_qubit_code()
    TN_code1 = TNCode(code)
    @test verify_code(TN_code1)
end

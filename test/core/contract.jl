using TensorNetworkCodes
using Test

@testset "contract" begin
    TN_code1 = TNCode(five_qubit_code())
    TN_code2 = TNCode(steane_code())
    TN_code12 = contract(TN_code1,TN_code2,[[1,2],[2,7]])
    @test verify_code(TN_code12)
end

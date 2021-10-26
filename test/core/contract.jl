using TensorNetworkCodes
using Test

@testset "contract" begin
    TN_code1 = TensorNetworkCode(five_qubit_code())
    TN_code2 = TensorNetworkCode(steane_code())
    TN_code12 = contract(TN_code1,TN_code2,[[1,2],[2,7]])
    @test verify_code(TN_code12)
end

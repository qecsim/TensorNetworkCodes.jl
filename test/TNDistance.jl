using TensorNetworkCodes
using TensorNetworkCodes.TNDistance
using Test

@testset "tn_distance" begin
    code = TensorNetworkCode(steane_code())
    code = contract(code, code, [[1, 1], [3, 3]])
    @test tn_distance(code) == 2
end

@testset "tn_operator_weights" begin
    code = TensorNetworkCode(five_qubit_code())
    code = contract(code, code, [[1, 1], [3, 3]])
    expected = OperatorWeights([1, 0, 0, 0, 9, 0, 6], [1, 0, 9, 24, 99, 72, 51])
    actual = tn_operator_weights(code)
    @test actual.stabilizer_weights == expected.stabilizer_weights
    @test actual.all_operator_weights == expected.all_operator_weights
    @test actual.distance == expected.distance
end

using TensorNetworkCodes
using TensorNetworkCodes.TNDistance
using Test

@testset "tn_distance" begin
    code = TensorNetworkCode(steane_code())
    code = contract(code, code, [[1, 1], [3, 3]])
    @test tn_distance(code) == 2
end

@testset "tn_distance-edge-case" begin
    # this case previously led to an error because after applying `contract`
    # one of the nodes has no "physical" indices
    tn_five_qubit = TensorNetworkCode(five_qubit_code())
    bell_pair = SimpleCode("Bell pair", [[1,1], [3,3]], [])
    tn_bell = TensorNetworkCode(bell_pair)
    code = contract(tn_five_qubit, tn_bell, [[1,1], [2,2]])
    @test tn_distance(code) == 1
end

@testset "tn_distance-contract-same-code" begin
    code = TensorNetworkCode(steane_code())
    code = contract(code, code, [[1, 1], [3, 3]])
    code = contract(code, code, [[1, 1], [7, 7]])
    brute_distance = find_distance_logicals(code)[1]
    tensor_distance = tn_distance(code)
    @test tensor_distance == brute_distance == 2
end

@testset "tn_operator_weights" begin
    code = TensorNetworkCode(five_qubit_code())
    code = contract(code, code, [[1, 1], [3, 3]])
    expected = OperatorWeights([1, 0, 0, 0, 9, 0, 6], [1, 0, 9, 24, 99, 72, 51])
    actual = tn_operator_weights(code)
    @test actual.stabilizer_weights == expected.stabilizer_weights
    @test actual.all_operator_weights == expected.all_operator_weights
    @test actual.distance == expected.distance
    plot_operator_weights(actual)  # no error thrown
end

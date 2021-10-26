using TensorNetworkCodes
using Test

@testset "SimpleCode" begin
    stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]]
    logicals = [[1, 1, 1, 1, 1], [3, 3, 3, 3, 3]]  # XXXXX, ZZZZZ
    pure_errors = [[0, 1, 0, 0, 0], [0, 0, 0, 0, 3], [0, 0, 3, 0, 0], [1, 0, 0, 0, 0]]
    code = SimpleCode("5-qubit code", stabilizers, logicals, pure_errors)
    @test verify_code(code)

    code = SimpleCode("5-qubit code", stabilizers, logicals)
    @test verify_code(code)

    code = SimpleCode()
    @test verify_code(code)
end

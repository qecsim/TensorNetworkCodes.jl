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

@testset "TensorNetworkCode" begin
    # construct TensorNetworkCode from SimpleCode
    simple_code = five_qubit_code()
    tn_code = TensorNetworkCode(simple_code)
    @test tn_code.stabilizers == simple_code.stabilizers
    @test tn_code.logicals == simple_code.logicals
    @test verify_code(tn_code)
    @test tn_code.seed_codes[simple_code.name] == simple_code
    # construct SimpleCode from TensorNetworkCode
    new_simple_code = SimpleCode(tn_code)
    @test new_simple_code.stabilizers == simple_code.stabilizers
    @test new_simple_code.logicals == simple_code.logicals
    @test verify_code(new_simple_code)
end

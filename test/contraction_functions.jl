using TensorNetworkCodes
using Test

@testset "combine" begin
    # simple codes
    code1 = five_qubit_code()
    code2 = steane_code()
    combined_code = combine(code1, code2)
    @test verify(combined_code)
    @test num_qubits(combined_code) == num_qubits(code1) + num_qubits(code2)
    @test length(combined_code.logicals) == length(code1.logicals) + length(code2.logicals)

    # tn codes
    code1 = TensorNetworkCode(five_qubit_code())
    code2 = TensorNetworkCode(steane_code())
    combined_code = combine(code1, code2)
    @test verify(combined_code)
    @test num_qubits(combined_code) == num_qubits(code1) + num_qubits(code2)
    @test length(combined_code.logicals) == length(code1.logicals) + length(code2.logicals)

    # compared combined_codes
    code1 = TensorNetworkCode(combine(five_qubit_code(), steane_code()))
    code2 = combine(TensorNetworkCode(five_qubit_code()), TensorNetworkCode(steane_code()))
    @test code1.stabilizers == code2.stabilizers
    @test code1.logicals == code2.logicals
end

@testset "fusion" begin
    # simple code
    code = combine(five_qubit_code(), steane_code())
    @test verify(code)
    # one-step fusion
    fused_code1 = fusion(code, [[1, 7], [2, 12]])
    @test verify(fused_code1)
    # step-by-step fusion
    fused_code2 = fusion(code, [[1, 7]])
    fused_code2 = fusion(fused_code2, [[1, 10]])
    @test verify(fused_code2)
    @test fused_code1.stabilizers == fused_code2.stabilizers
    @test num_qubits(fused_code1) == num_qubits(fused_code1) == num_qubits(code) - 4

    # tn code
    code = combine(TensorNetworkCode(five_qubit_code()), TensorNetworkCode(steane_code()))
    @test verify(code)
    # one-step fusion
    fused_code1 = fusion(code, [[1, 7], [2, 12]])
    @test verify(fused_code1)
    # step-by-step fusion
    fused_code2 = fusion(code, [[1, 7]])
    fused_code2 = fusion(fused_code2, [[1, 10]])
    @test verify(fused_code2)
    @test fused_code1.stabilizers == fused_code2.stabilizers
    @test num_qubits(fused_code1) == num_qubits(fused_code1) == num_qubits(code) - 4

    # STRANGE MESSAGE: "Self contraction occurred.  Contraction algorithms may not work!"
    # NOTE: swap the "code = ..." lines below to see the message.
    # code = TensorNetworkCode(combine(five_qubit_code(), steane_code()))
    code = combine(TensorNetworkCode(five_qubit_code()), TensorNetworkCode(steane_code()))
    fusion(code, [[1, 7], [2, 12]])
end

@testset "contract" begin
    # tn code
    code1 = TensorNetworkCode(five_qubit_code())
    code2 = TensorNetworkCode(steane_code())
    contracted_code = contract(code1, code2, [[1, 2], [2, 7]])
    @test verify(contracted_code)
    # equivalent using combine and fusion
    equivalent_code = fusion(combine(code1, code2), [[1, 7], [2, 12]])
    @test contracted_code.stabilizers == equivalent_code.stabilizers
    @test contracted_code.logicals == equivalent_code.logicals
end

@testset "contract_by_coords" begin
    code1 = TensorNetworkCode(five_qubit_code())
    code2 = TensorNetworkCode(steane_code())
    # set coordinates so contract on [1, 2] and [2, 7] according to contract method
    code1_coords = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5]]
    code2_coords = [[1, 0], [1, 1], [0, 1], [1, 3], [1, 4], [1, 5], [1, 6], [0, 2]]
    set_coords!(code1, code1_coords)
    set_coords!(code2, code2_coords)
    contracted_code = contract_by_coords(code1, code2)
    @test verify(contracted_code)
    # equivalent using contract
    equivalent_code = contract(code1, code2, [[1, 2], [2, 7]])
    @test contracted_code.stabilizers == equivalent_code.stabilizers
    @test contracted_code.logicals == equivalent_code.logicals
end

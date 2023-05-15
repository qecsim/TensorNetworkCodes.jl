using TensorNetworkCodes
using TensorNetworkCodes.TNDecode
using Test

@testset "tn_decode" begin
    # surface code
    code = rotated_surface_code(3)
    p = 0.2;  # error_probability
    error = pauli_random_operator(num_qubits(code), p)
    syndrome = find_syndrome(code, error);
    # default_contract
    recovery, success_prob = tn_decode(code, syndrome, p)
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1
    # basic_contract (same as default)
    recovery, success_prob = tn_decode(code, syndrome, p; contract_fn=basic_contract())
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1
    # mps_contract
    recovery, success_prob = tn_decode(code, syndrome, p; contract_fn=mps_contract(8))
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1

    # small non-uniform code
    tn_five_qubit = TensorNetworkCode(five_qubit_code())
    tn_six_qubit = TensorNetworkCode(purify(five_qubit_code()))
    bell_pair = SimpleCode("Bell pair", [[1,1], [3,3]], [])
    tn_bell = TensorNetworkCode(bell_pair)
    code = TensorNetworkCodes.contract(tn_five_qubit, tn_six_qubit, [[1,1],[2,2]])
    code = TensorNetworkCodes.contract(code, tn_bell, [[1,1], [2,2]])
    p = 0.2;  # error_probability
    error = pauli_random_operator(num_qubits(code), p)
    syndrome = find_syndrome(code, error)
    # default_contract
    recovery, success_prob = tn_decode(code, syndrome, p)
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1
    # basic_contract (same as default)
    recovery, success_prob = tn_decode(code, syndrome, p; contract_fn=basic_contract())
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1
    # mps_contract
    recovery, success_prob = tn_decode(code, syndrome, p; contract_fn=mps_contract(8))
    @test find_syndrome(code, recovery) == syndrome
    @test 0 <= success_prob <= 1
end

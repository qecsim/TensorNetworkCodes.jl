using TensorNetworkCodes
using TensorNetworkCodes.TNDecode
using Test

@testset "tn_decode" begin
    code = rotated_surface_code(3)
    p = 0.2;  # error_probability
    error = random_pauli_error(num_qubits(code), p)
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
end

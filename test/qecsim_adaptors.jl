using Qecsim
using Qecsim.GenericModels: DepolarizingErrorModel, NaiveDecoder
using TensorNetworkCodes
using TensorNetworkCodes: _bsf_to_tnpauli, _tnpauli_to_bsf
using Test

@testset "_tnpauli_to_bsf" begin
    # single operator
    str_pauli = "IXIYIZ"
    int_pauli = pauli_rep_change.(collect(str_pauli))
    bsf_pauli = _tnpauli_to_bsf(int_pauli)
    @test bsf_pauli == Qecsim.to_bsf(str_pauli)
    # multiple operators
    str_paulis = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    int_paulis = [pauli_rep_change.(collect(p)) for p in str_paulis]
    bsf_paulis = _tnpauli_to_bsf(int_paulis)
    @test bsf_paulis == Qecsim.to_bsf(str_paulis)
end

@testset "_bsf_to_tnpauli" begin
    # single operator
    str_pauli = "IXIYIZ"
    bsf_pauli = Qecsim.to_bsf(str_pauli)
    int_pauli = _bsf_to_tnpauli(bsf_pauli)
    @test int_pauli == pauli_rep_change.(collect(str_pauli))
    # multiple operators
    str_paulis = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    bsf_paulis = Qecsim.to_bsf(str_paulis)
    int_paulis = _bsf_to_tnpauli(bsf_paulis)
    @test int_paulis == [pauli_rep_change.(collect(p)) for p in str_paulis]
end

@testset "QecsimTNCode" begin
    # default distance, label
    tn_code = TensorNetworkCode(five_qubit_code())
    qs_code = QecsimTNCode(tn_code)
    @test validate(qs_code) === nothing  # no error
    @test isequal(nkd(qs_code), (5, 1, missing))
    @test label(qs_code) == "QecsimTNCode: [5,1,missing]"
    # kwargs distance, label
    tn_code = TensorNetworkCode(steane_code())
    qs_code = QecsimTNCode(tn_code; distance=3, label="Steane")
    @test validate(qs_code) === nothing  # no error
    @test isequal(nkd(qs_code), (7, 1, 3))
    @test label(qs_code) == "Steane"
end

@testset "QecsimTNDecoder" begin
    # field test
    @test label(QecsimTNDecoder()) == "QecsimTNDecoder"
    @test label(QecsimTNDecoder(1)) == "QecsimTNDecoder (chi=1)"
    # invalid parameters
    @test_throws ArgumentError QecsimTNDecoder(-1)
    # models
    tn_code = rotated_surface_code(3)
    qs_code = QecsimTNCode(tn_code; distance=3, label="Rotated Surface 3")
    qs_error_model = DepolarizingErrorModel()
    qs_decoder = QecsimTNDecoder()
    p = 0.1
    # direct test of decoder (exact contraction)
    qs_error = generate(qs_error_model, qs_code, p)
    qs_syndrome = bsp(stabilizers(qs_code), qs_error)
    qs_result = decode(qs_decoder, qs_code, qs_syndrome; p=p)
    @test bsp(stabilizers(qs_code), qs_result.recovery) == qs_syndrome
    @test !any(bsp(stabilizers(qs_code), xor.(qs_error, qs_result.recovery)))
    @test 0 <= qs_result.custom_values[1][1] <= 1 # success probability
    # direct test of decoder (approx. contraction)
    qs_result = decode(QecsimTNDecoder(1), qs_code, qs_syndrome; p=p)
    @test bsp(stabilizers(qs_code), qs_result.recovery) == qs_syndrome
    @test !any(bsp(stabilizers(qs_code), xor.(qs_error, qs_result.recovery)))
    @test 0 <= qs_result.custom_values[1][1] <= 1 # success probability
    # test via run_once
    qec_run_once(qs_code, qs_error_model, qs_decoder, p)  # no error
end

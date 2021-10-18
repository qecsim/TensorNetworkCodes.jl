using TensorNetworkCodes
using Test

@testset "pauli_commutation" begin
    # single-qubit
    @test pauli_commutation.(Ref(0), [0, 1, 2, 3]) == [0, 0, 0, 0]  # I with IXYZ
    @test pauli_commutation.(Ref(1), [0, 1, 2, 3]) == [0, 0, 1, 1]  # X with IXYZ
    @test pauli_commutation.(Ref(2), [0, 1, 2, 3]) == [0, 1, 0, 1]  # Y with IXYZ
    @test pauli_commutation.(Ref(3), [0, 1, 2, 3]) == [0, 1, 1, 0]  # Z with IXYZ
    # multi-qubit
    # 5-qubit code stabilizers: XZZXI, IXZZX, XIXZZ, ZXIXZ
    stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]]
    @test pauli_commutation.(Ref([1, 1, 1, 1, 1]), stabilizers) == zeros(4) # XXXXX
    @test pauli_commutation.(Ref([3, 3, 3, 3, 3]), stabilizers) == zeros(4) # ZZZZZ
    @test pauli_commutation.(Ref([1, 1, 0, 0, 0]), stabilizers) == [1, 0, 0, 1] # XXIII
    @test pauli_commutation.(Ref([0, 0, 0, 3, 3]), stabilizers) == [1, 1, 0, 1] # IIIZZ
    @test pauli_commutation.(Ref([0, 0, 2, 0, 0]), stabilizers) == [1, 1, 1, 0] # IIYII
end

@testset "pauli_product" begin
    @test pauli_product.(Ref(0), [0, 1, 2, 3]) == [0, 1, 2, 3]  # I with IXYZ -> IXYZ
    @test pauli_product.(Ref(1), [0, 1, 2, 3]) == [1, 0, 3, 2]  # X with IXYZ -> XIZY
    @test pauli_product.(Ref(2), [0, 1, 2, 3]) == [2, 3, 0, 1]  # Y with IXYZ -> YZIX
    @test pauli_product.(Ref(3), [0, 1, 2, 3]) == [3, 2, 1, 0]  # Z with IXYZ -> ZYXI
end

@testset "are_they_independent" begin
    new_code = steane_code()
    @test are_they_independent(new_code.stabilizers)
end

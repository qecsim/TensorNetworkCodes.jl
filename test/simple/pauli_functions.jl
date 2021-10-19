using TensorNetworkCodes
using Test

@testset "pauli_are_commuting" begin
    # 5-qubit code stabilizers: XZZXI, IXZZX, XIXZZ, ZXIXZ
    stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]]
    logical_x = [1, 1, 1, 1, 1]  # XXXXX
    logical_z = [3, 3, 3, 3, 3]  # ZZZZZ
    @test pauli_are_commuting(stabilizers)
    @test pauli_are_commuting(vcat(stabilizers, [logical_x]))
    @test pauli_are_commuting(vcat(stabilizers, [logical_z]))
    @test !pauli_are_commuting([logical_x, logical_z])
end

@testset "pauli_are_independent" begin
    # 5-qubit code stabilizers: XZZXI, IXZZX, XIXZZ, ZXIXZ
    stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]]
    logical_x = [1, 1, 1, 1, 1]  # XXXXX
    logical_z = [3, 3, 3, 3, 3]  # ZZZZZ
    @test pauli_are_independent(stabilizers)
    @test pauli_are_independent(vcat(stabilizers, [logical_x, logical_z]))
    extra_stabilizer = pauli_product.(stabilizers[1], stabilizers[2])
    @test !pauli_are_independent(vcat(stabilizers, [extra_stabilizer]))
end

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

@testset "pauli_product_pow" begin
    ops = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1]]
    @test pauli_product_pow(ops, [0, 0, 0]) == [0, 0, 0, 0]
    @test pauli_product_pow(ops, [0, 0, 1]) == [2, 3, 0, 1]
    @test pauli_product_pow(ops, [0, 1, 0]) == [1, 0, 3, 2]
    @test pauli_product_pow(ops, [0, 1, 1]) == [3, 3, 3, 3]
    @test pauli_product_pow(ops, [1, 0, 0]) == [0, 1, 2, 3]
    @test pauli_product_pow(ops, [1, 0, 1]) == [2, 2, 2, 2]
    @test pauli_product_pow(ops, [1, 1, 0]) == [1, 1, 1, 1]
    @test pauli_product_pow(ops, [1, 1, 1]) == [3, 2, 1, 0]
end

@testset "pauli_rep_change" begin
    @test pauli_rep_change.([0, 1, 2, 3]) == ['I', 'X', 'Y', 'Z']
    @test pauli_rep_change.(['I', 'X', 'Y', 'Z']) == [0, 1, 2, 3]
end

@testset "pauli_weight" begin
    @test pauli_weight([0, 0, 0, 0, 0]) == 0
    @test pauli_weight([0, 0, 1, 2, 0]) == 2
    @test pauli_weight([3, 0, 2, 0, 1]) == 3
    @test pauli_weight([3, 2, 1, 2, 2]) == 5
end

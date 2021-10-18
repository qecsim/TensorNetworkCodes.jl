using TensorNetworkCodes
using Test

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

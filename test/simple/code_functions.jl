using TensorNetworkCodes
using Test

@testset "num_qubits" begin
    @test num_qubits(five_qubit_code()) == 5
    @test num_qubits(steane_code()) == 7
    @test num_qubits(random_code(4, 1)) == 4
    @test num_qubits(random_stabilizer_state(6)) == 6
end

@testset "verify_code" begin
    @test verify_code(random_code(6,2))
end

@testset "gauge_code" begin
    new_code = five_qubit_surface_code()
    new_code = gauge_code(new_code,[[0,1]],[1])
    @test verify_code(new_code)
end

using TensorNetworkCodes
using Test

@testset "five_qubit_code" begin
    code = five_qubit_code()
    @test verify(code)
    @test num_qubits(code) == 5
    @test Int(length(code.logicals) / 2) == 1
    @test find_distance_logicals(code)[1] == 3
end

@testset "five_qubit_surface_code" begin
    code = five_qubit_surface_code()
    @test verify(code)
    @test num_qubits(code) == 5
    @test Int(length(code.logicals) / 2) == 1
    @test find_distance_logicals(code)[1] == 2
end

@testset "steane_code" begin
    code = steane_code()
    @test verify(code)
    @test num_qubits(code) == 7
    @test Int(length(code.logicals) / 2) == 1
    @test find_distance_logicals(code)[1] == 3
end

@testset "random_code" begin
    code = random_code(9, 2)
    @test verify(code)
    @test num_qubits(code) == 9
    @test Int(length(code.logicals) / 2) == 2
    # invalid k > n
    @test_throws ErrorException random_code(5, 6)
end

@testset "random_stabilizer_state" begin
    code = random_stabilizer_state(6)
    @test verify(code)
    @test num_qubits(code) == 6
    @test Int(length(code.logicals) / 2) == 0
end

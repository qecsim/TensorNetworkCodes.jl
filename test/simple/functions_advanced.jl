using TensorNetworkCodes
using Test

@testset "are_they_independent" begin
    new_code = steane_code()
    @test are_they_independent(new_code.stabilizers)
end

@testset "verify_code" begin
    @test verify_code(random_code(6,2))
end

@testset "gauge_code" begin
    new_code = five_qubit_surface_code()
    new_code = gauge_code(new_code,[[0,1]],[1])
    @test verify_code(new_code)
end

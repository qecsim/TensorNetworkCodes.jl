using SimpleCodes
using Test

@testset "SimpleCodes.jl" begin
    @test size(five_qubit_code()) == 5
end


@testset "Code_functions_advanced.jl" begin
    new_code = steane_code()
    @test are_they_independent(new_code.stabilizers)
    @test verify_code(random_code(6,2))

    new_code = five_qubit_surface_code()
    new_code = gauge_code(new_code,[[0,1]],[1])
    @test verify_code(new_code)
end
using TensorNetworkCodes
using Test

@testset "size" begin
    @test size(five_qubit_code()) == 5
end

@testset "verify_code" begin
    @test verify_code(random_code(6,2))
end

@testset "gauge_code" begin
    new_code = five_qubit_surface_code()
    new_code = gauge_code(new_code,[[0,1]],[1])
    @test verify_code(new_code)
end

using TensorNetworkCodes
using Test

@testset "size" begin
    @test size(five_qubit_code()) == 5
end

using TensorNetworkCodes
using Test

@testset "are_they_independent" begin
    new_code = steane_code()
    @test are_they_independent(new_code.stabilizers)
end

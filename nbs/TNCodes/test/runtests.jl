using SimpleCodes
using TNCodes
using Test

@testset "TNCodes.jl" begin
    # test TN_code_types.jl
    code = five_qubit_code()
    TN_code1 = TNCode(code)
    @test verify_code(TN_code1)

    # test TN_code_contract.jl
    TN_code2 = TNCode(steane_code())
    TN_code12 = contract(TN_code1,TN_code2,[[1,2],[2,7]])
    @test verify_code(TN_code12)

    # test TN_surface.jl
    surface = surface_code(3)
    @test verify_code(surface)
    @test distance(surface)[1] == 4

end
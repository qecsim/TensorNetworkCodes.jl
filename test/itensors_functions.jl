using TensorNetworkCodes
using TensorNetworkCodes: _code_to_tensor
using Test
using ITensors

@testset "_code_to_tensor" begin
    # Test on the five-qubit code
    tensor = _code_to_tensor(five_qubit_code())
    # Check total number of code operators is correct
    @test sum(tensor) == 64.0
    # Check dimensions are correct
    @test size(tensor) == (4, 4, 4, 4, 4, 4)

    #Repeat for Steane code
    tensor = _code_to_tensor(steane_code())
    @test sum(tensor) == 256.0
    @test size(tensor) == (4, 4, 4, 4, 4, 4, 4, 4)
end

@testset "code_to_Itensor" begin
    logical_indices = [Index(4,"logical")]
    physical_indices = [Index(4,"physical") for _ in 1:5]
    tensor = code_to_Itensor(five_qubit_surface_code(),logical_indices,physical_indices)
    # Check dimensions
    @test dims(tensor) == (4, 4, 4, 4, 4, 4)
end

@testset "identity_coset" begin
    logical_indices = [Index(4,"logical")]
    physical_indices = [Index(4,"physical") for _ in 1:7]
    tensor = code_to_Itensor(steane_code(),logical_indices,physical_indices)
    id_coset = identity_coset(tensor)
    # Check number of stabilizers is correct
    @test sum(id_coset) == 64.0
end

@testset "all_cosets" begin
    logical_indices = [Index(4,"logical")]
    physical_indices = [Index(4,"physical") for _ in 1:5]
    tensor = code_to_Itensor(five_qubit_surface_code(),logical_indices,physical_indices)
    all_tensor = all_cosets(tensor)
    # Check number of stabilizers + logicals is correct
    @test sum(all_tensor) == 64.0
end

@testset "physical_tensor" begin
    index = Index(4,"logical")
    error_prob = 0.3
    pauli = 3
    tensor = physical_tensor(index,error_prob,pauli)
    arr = array(tensor)
    @test round.(arr,digits=1) == [0.1, 0.1, 0.1, 0.7]
end

@testset "create_virtual_tensor" begin
    code = five_qubit_surface_code()
    tn_code = TensorNetworkCode(code)
    tn_code = TensorNetworkCodes.contract(tn_code,deepcopy(tn_code),[[1,1]])
    tensor = create_virtual_tensor(tn_code,-1)
    # Check to see if tensor has all the right kind of indices
    @test hastags(tensor,"physical")
    @test hastags(tensor,"logical")
    @test hastags(tensor,"bond")
end
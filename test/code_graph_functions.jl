using TensorNetworkCodes
using Test

@testset "edges" begin
    code = TensorNetworkCode(five_qubit_code())
    @test Set(edges(code)) == Set([Set([-1, i]) for i in 1:5])
end

@testset "nodes" begin
    code = TensorNetworkCode(five_qubit_code())
    @test nodes(code) == [-1, 1, 2, 3, 4, 5]
end

@testset "num_nodes" begin
    code = TensorNetworkCode(five_qubit_code())
    @test num_nodes(code) == 6  # 5 physical and 1 logical qubits
end

@testset "physical_neighbours" begin
    code = TensorNetworkCode(five_qubit_code())
    @test physical_neighbours(code, -1) == Set(1:5)
    @test physical_neighbours(code, 2) == Set()
end

@testset "set_coords!" begin
    code = TensorNetworkCode(five_qubit_code())
    old_coords = deepcopy(coords(code))
    @test coords(code) == old_coords
    # let's put them in a line
    new_coords = [[0, i] for i in 1:num_nodes(code)]
    set_coords!(code, new_coords)
    @test coords(code) != old_coords
    @test [coords(code, n) for n in nodes(code)] == new_coords
end

@testset "shift_coords!" begin
    code = TensorNetworkCode(five_qubit_code())
    old_coords = deepcopy(coords(code))
    @test coords(code) == old_coords
    # let's shift a bit
    shift = [1, 2]
    shift_coords!(code, shift)
    @test coords(code) != old_coords
    @test coords(code) == Dict([k, v + shift] for (k, v) in old_coords)
end

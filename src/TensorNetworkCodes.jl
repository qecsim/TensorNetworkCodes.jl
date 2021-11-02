module TensorNetworkCodes

# imports
using Combinatorics: combinations
using GraphPlot: gplot
using ITensors: Index, ITensor, dim, hastags, inds, tags
using LightGraphs: Graph, add_edge!
using Random: AbstractRNG, MersenneTwister, RandomDevice, GLOBAL_RNG, rand
using StatsBase: ProbabilityWeights, sample

# exports
export QuantumCode, SimpleCode, CodeGraph, TensorNetworkCode
include("types.jl")
export pauli_are_commuting, pauli_are_independent, pauli_commutation
export pauli_pow, pauli_product, pauli_product_pow
export pauli_random_operator, pauli_rep_change, pauli_weight
include("pauli_functions.jl")
export num_qubits, verify
export find_distance_logicals, find_pure_error, find_pure_errors, find_syndrome
export gauge, permute, purify
include("code_functions.jl")
export edges, nodes, num_nodes, physical_neighbours, set_coords!, shift_coords!
export coords, edge_indices, edge_types, node_indices, node_types
export set_coords!, set_edge_indices!, set_edge_types!, set_node_indices!, set_node_types!, new_indices
include("code_graph_functions.jl")
export combine, contract, contract_by_coords, fusion
include("contraction_functions.jl")
export plot_code, plot_operator
include("plotting_functions.jl")
export five_qubit_code, five_qubit_surface_code, steane_code
export random_code, random_stabilizer_state
include("examples/simple.jl")
export almost_rotated_surface_code, rotated_surface_code, surface_code
include("examples/surface.jl")
export code_to_Itensor, identity_coset, all_cosets
export physical_tensor, create_virtual_tensor
include("itensors_functions.jl")

# include submodules (not exported)
include("TNDecode.jl")
include("TNDistance.jl")
include("QecsimAdaptors.jl")

end

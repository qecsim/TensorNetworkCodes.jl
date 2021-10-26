module TensorNetworkCodes

# imports
using Qecsim

# SimpleCode imports
using Combinatorics: combinations
using Random: MersenneTwister, RandomDevice, rand
using StatsBase: Weights, sample

# Core imports
using GraphPlot: gplot
using ITensors: Index, ITensor, dim, hastags, inds
using LightGraphs: Graph, add_edge!

# exports
export QuantumCode, SimpleCode, CodeGraph, TensorNetworkCode
include("types.jl")
export pauli_are_commuting, pauli_are_independent, pauli_commutation
export pauli_pow, pauli_product, pauli_product_pow, pauli_rep_change, pauli_weight
include("pauli_functions.jl")
export num_qubits, verify_code
export find_distance_logicals, find_pure_error, find_pure_errors, find_syndrome
export gauge_code, permute_code, purify_code
include("code_functions.jl")
export edges, nodes, num_nodes, shift_coords!
export coords, edge_indices, edge_types, node_indices, node_types
export set_coords!, set_edge_indices!, set_edge_types!, set_node_indices!, set_node_types!
include("code_graph_functions.jl")

# SimpleCode exports
export five_qubit_code,five_qubit_surface_code,steane_code
export random_code, random_stabilizer_state
include("simple/examples.jl")
export random_pauli_error
include("simple/errors.jl")
export min_weight_brute_force,do_nothing_decoder
export monte_carlo_simulation
include("simple/decoding.jl")

# Core exports
export identity_coset,all_cosets,gauge_code,code_to_tensor,code_to_Itensor
include("core/functions.jl")
export code_plot,operator_plot
include("core/plotting.jl")
export combine_by_coordinates,contract,fusion,merge
include("core/contract.jl")
export surface_code,rotated_surface_code,almost_rotated_surface_code
include("core/surface.jl")

# include submodules (not reexported)
include("TNDecode.jl")
include("TNDistance.jl")

end

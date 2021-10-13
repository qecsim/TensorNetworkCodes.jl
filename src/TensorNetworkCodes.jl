module TensorNetworkCodes

# imports
using Qecsim

# SimpleCode imports
import Base: size
import StatsBase
using Combinatorics
using Random

# Core imports
import Base.merge
using LightGraphs, GraphPlot, ITensors


# exports

# SimpleCode exports
export SimpleCode, QuantumCode
include("simple/types.jl")
export five_qubit_code,five_qubit_surface_code,steane_code
include("simple/examples.jl")
export size, pauli_product, do_they_commute
export pauli_rep_change, permute, weight
include("simple/functions.jl")
export are_they_independent, generate_pure_errors, verify_code
export distance, low_weight_stabilizers, are_physically_equivalent
export purify_code, gauge_code, random_stabilizer_state, random_code
include("simple/functions_advanced.jl")
export random_pauli_error
include("simple/errors.jl")
export get_syndrome,get_pure_error,min_weight_brute_force,do_nothing_decoder
export monte_carlo_simulation
include("simple/decoding.jl")

# Core exports
export TNCode,num_nodes,nodes,edges,SimpleCode
export coords,node_types,node_indices,edge_types,edge_indices
export set_coords!,set_node_types!,set_node_indices!,set_edge_types!,set_edge_indices!
export shift_coords!
include("core/types.jl")
export identity_coset,all_cosets,gauge_code,code_to_tensor,code_to_Itensor
include("core/functions.jl")
export code_plot,operator_plot
include("core/plotting.jl")
export combine_by_coordinates,contract,fusion,merge
include("core/contract.jl")
export surface_code,diamond_lattice_code,almost_surface_code,fully_random_code
export checkerboard_code, funny_code
include("core/surface.jl")

end

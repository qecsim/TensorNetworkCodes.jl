module TNCodes

# imports
using SimpleCodes
import SimpleCodes.gauge_code


# exports
export TNCode,num_nodes,nodes,edges,SimpleCode
export coords,node_types,node_indices,edge_types,edge_indices
export set_coords!,set_node_types!,set_node_indices!,set_edge_types!,set_edge_indices!
export shift_coords!
include("TN_code_types.jl")

export identity_coset,all_cosets,gauge_code,code_to_tensor,code_to_Itensor
include("TN_code_functions.jl")

export code_plot,operator_plot
include("TN_code_plotting.jl")

export combine_by_coordinates,contract,fusion,merge
include("TN_code_contract.jl")

export surface_code,diamond_lattice_code,almost_surface_code,fully_random_code
export checkerboard_code, funny_code
include("TN_surface.jl")

end

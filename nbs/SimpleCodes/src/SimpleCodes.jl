module SimpleCodes


# imports

# exports
export SimpleCode, QuantumCode
include("Code_types.jl")

export five_qubit_code,five_qubit_surface_code,steane_code
include("Code_examples.jl")

export size, pauli_product, do_they_commute
export pauli_rep_change, permute, weight
include("Code_functions.jl")

export are_they_independent, generate_pure_errors, verify_code
export distance, low_weight_stabilizers, are_physically_equivalent
export purify_code, gauge_code, random_stabilizer_state, random_code
include("Code_functions_advanced.jl")

export random_pauli_error
include("Code_errors.jl")

export get_syndrome,get_pure_error,min_weight_brute_force,do_nothing_decoder
export monte_carlo_simulation
include("Code_decoding.jl")

end

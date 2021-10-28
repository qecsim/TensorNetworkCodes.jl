"""
    code_to_tensor(code)

Returns a tensor (array) describing the logical cosets of the code.
"""
function code_to_tensor(code::QuantumCode)

    alt_code = purify(code)
    n = num_qubits(alt_code)
    g = alt_code.stabilizers

    dims = fill(4,n)
    tensor = zeros(Float64,dims...)

    for α in 0:2^n-1
        powers = digits!(zeros(Int64,n),α,base = 2)
        operator = pauli_product_pow(g,powers)
        operator = operator .+ 1

        index = CartesianIndex(operator...)
        tensor[index] = 1
    end

    return tensor
end





"""
    code_to_Itensor(code,indices)

Returns an `ITensor` with indices describing the logical cosets of the code.
"""
function code_to_Itensor(code::QuantumCode,
        logical_indices::Array{Index{Int64},1},
        physical_indices::Array{Index{Int64},1})

    n = num_qubits(code)
    K = length(code.logicals)
    if K/2 != length(logical_indices) || n != length(physical_indices)
        error("incorrect number of indices!")
    end

    tensor = code_to_tensor(code)
    indices = vcat(logical_indices,physical_indices)

    return ITensor(tensor,indices...)
end





"""
    identity_coset(tensor)

Given an `ITensor` descibing a code, return the `ITensor` describing
only the identity coset, i.e., the stabilizer group.
"""
function identity_coset(tensor::ITensor)
    indices = inds(tensor)
    logical_indices = [χ for χ in indices if hastags(χ,"logical")]

    for index in logical_indices
        ρ = ITensor(index)
        ρ[index => 1] = 1
        tensor = tensor * ρ
    end

    return tensor
end





"""
    all_cosets(tensor)

Given an `ITensor` descibing a code, return the `ITensor` describing
all the cosets by summing over all logical indices.
"""
function all_cosets(tensor::ITensor)
    indices = inds(tensor)
    logical_indices = [χ for χ in indices if hastags(χ,"logical")]

    for index in logical_indices
        ρ = ITensor(ones(dim(index)),index)
        tensor = tensor * ρ
    end

    return tensor
end

"""
    physical_tensor(index,error_prob,pauli)

Returns a one leg `ITensor` describing the error on one site.
"""
function physical_tensor(
        index::Index{Int64},
        error_prob::Float64,
        pauli::Int64)

    array = fill(error_prob/3,4)
    array[pauli + 1] = 1 - error_prob

    return ITensor(array,index)
end

"""
    create_virtual_tensor(graph,node)

Creates the `ITensor` describing the seed code at `node`.
"""
function create_virtual_tensor(code::TensorNetworkCode,node::Int64)

    indices = node_indices(code,node)
    code_type = node_types(code,node)
    seed_code = code.seed_codes[code_type]

    k = num_qubits(seed_code) - length(seed_code.stabilizers)
    logical_indices = [Index(4,"logical") for α in 1:k]

    return code_to_Itensor(
        seed_code,
        logical_indices,
        indices)
end
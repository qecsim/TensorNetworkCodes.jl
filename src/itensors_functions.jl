"""
    _code_to_tensor(code::QuantumCode) -> Array{Float64,N} where N

Given a [`QuantumCode`](@ref), returns a tensor describing the 
logical cosets of the code.
"""
function _code_to_tensor(code::QuantumCode)

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
    code_to_Itensor(code::QuantumCode,logical_indices::Array{Index{Int64},1},physical_indices::Array{Index{Int64},1})
    -> ITensor

Returns an `ITensor` describing the logical cosets of the code.  Indices 
(`Index{Int64}`) correspond to logical and physical qubits of the code.

# Examples
```jldoctest
julia> using ITensors;

julia> logical_indices = [Index(4,"logical")];

julia> physical_indices = [Index(4,"physical") for _ in 1:5];

julia> tensor = code_to_Itensor(five_qubit_surface_code(),logical_indices,physical_indices);

julia> dims(tensor) # tensor has six legs
(4, 4, 4, 4, 4, 4)
```
"""
function code_to_Itensor(code::QuantumCode,
        logical_indices::Array{Index{Int64},1},
        physical_indices::Array{Index{Int64},1})

    n = num_qubits(code)
    K = length(code.logicals)
    if K/2 != length(logical_indices) || n != length(physical_indices)
        error("incorrect number of indices!")
    end

    tensor = _code_to_tensor(code)
    indices = vcat(logical_indices,physical_indices)

    return ITensor(tensor,indices...)
end





"""
    identity_coset(tensor::ITensor) -> ITensor

Given an `ITensor` descibing a code, return the `ITensor` describing
only the identity coset, i.e., the stabilizer group.

# Examples
```jldoctest
julia> using ITensors;

julia> logical_indices = [Index(4,"logical")];

julia> physical_indices = [Index(4,"physical") for _ in 1:7];

julia> tensor = code_to_Itensor(steane_code(),logical_indices,physical_indices);

julia> sum(identity_coset(tensor)) # stabilizer group has 64 elements
64.0
```
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
    all_cosets(tensor::ITensor) -> ITensor

Given an `ITensor` descibing a code, return the `ITensor` describing
all the cosets by summing over all logical indices.  Works similarly to 
[`identity_coset`](@ref).
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
    physical_tensor(index::Index{Int64},error_prob::Float64,pauli::Int64) -> ITensor

Returns a one leg `ITensor` describing the error on one physical qubit.
Assumes depolarizing noise as the error model.  `pauli` permutes the values
(useful for making decoders to account for pure_errors).

# Examples
```jldoctest
julia> using ITensors;

julia> index = Index(4,"physical");

julia> pauli = 0;

julia> error_prob = 0.3;

julia> a = array(physical_tensor(index,error_prob,pauli));

julia> println(round.(a,digits=1))
[0.7, 0.1, 0.1, 0.1]

julia> pauli = 2;

julia> a = array(physical_tensor(index,error_prob,pauli));

julia> println(round.(a,digits=1))
[0.1, 0.1, 0.7, 0.1]
```
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
    create_virtual_tensor(code::TensorNetworkCode,node::Int64) -> ITensor

Creates the `ITensor` describing the seed code at `node`.  Necessary for, e.g.,
tensor-network decoding or distance calculation.
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
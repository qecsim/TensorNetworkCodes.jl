"""
Decoding functions based on contracting a tensor-network associated with the tensor-network
code, see [`tn_decode`](@ref).
"""
module TNDecode

#imports
using ..TensorNetworkCodes: TensorNetworkCode
using ..TensorNetworkCodes: edges, edge_indices, node_indices, node_types
using ..TensorNetworkCodes: code_to_Itensor, create_virtual_tensor, physical_tensor
using ..TensorNetworkCodes: pauli_commutation, pauli_product, pauli_product_pow
using ..TensorNetworkCodes: num_qubits, physical_neighbours
using ..TensorNetworkCodes: find_pure_error, find_syndrome
using ..TensorNetworkCodes: random_pauli_error
using ITensors: ITensors  # only imported to avoid `contract` name clash
using ITensors: Index, ITensor, MPO, MPS, array, hastags
using Statistics: mean, stdm

#exports
export basic_contract, mps_contract, tn_decode

# ########## Should be in core/types.jl
# function _physical_neighbours(code,node)

#     all_edges = edges(code)
#     neighbour_edges = filter(x->node in x,all_edges)
#     neighbour_nodes = union(neighbour_edges...)

#     return filter(x->x>0,neighbour_nodes)

# end
# ##########

"""
    basic_contract(bond_dim=10) -> function

Return an exact contract function for use with [`tn_decode`](@ref) that contracts a matrix
of `ITensor` working along row-by-row.

See also: [`tn_decode`](@ref).
"""
function basic_contract()
    return function _basic_contract(imatrix)
        output = ITensor(1)
        for i in 1:size(imatrix)[1], j in 1:size(imatrix)[2]
            output = output * imatrix[i,j]
        end

        return output
    end
end

"""
    mps_contract(bond_dim=10) -> function

Return a contract function for use with [`tn_decode`](@ref) that contracts a matrix of
`ITensor`, treating rows as MPS/MPO and contracting rows element-wise, retaining a maximum
bond dimension of `bond_dim`, before contracting along the final row.

See also: [`tn_decode`](@ref).
"""
function mps_contract(bond_dim=10)
    # let trick for captured variable performance
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
    fn = let bond_dim = bond_dim  #
        function _mps_contract(imatrix)
            ψ = MPS(imatrix[1,:])
            for j in 2:size(imatrix)[1]-1
                Σ = MPO(imatrix[j,:])
                ψ = ITensors.contract(Σ, ψ; maxdim = bond_dim)  # ITensors.contract
            end

            output = imatrix[end,1] * ψ[1]
            for j in 2:size(imatrix)[2]
                output = output * imatrix[end,j] * ψ[j]
            end

            return output
        end
    end
    return fn
end

"""
    tn_decode(
        code::TensorNetworkCode, syndrome::AbstractVector{Int}, error_probability::Real;
        contract_fn=basic_contract()
    ) -> (Vector{Int}, Float64)

Return a recovery operator for the code, that is consistent with the syndrome and is from
the most-likely logical coset. Additionally, the predicted success probability, defined as
the ratio of the probability of the mostly-likely logical coset to the sum of the
probabilities of all logical cosets, is returned.

The method of tensor-network contraction defaults to `basic_contract`, which corresponds to
exact contraction. Alternatively, `mps_contract` may be specified for efficient but
approximate contraction.

See also: [`basic_contract`](@ref), [`mps_contract`](@ref)

!!! note
    This function currently assumes a depolarizing noise model and only supports codes that
    can be laid out on a square lattice.

# Examples
```jldoctest
julia> using TensorNetworkCodes.TNDecode

julia> code = rotated_surface_code(3);

julia> error = [0, 0, 0, 1, 0, 0, 3, 0, 0];  # IIIXIIZII

julia> syndrome = find_syndrome(code, error);

julia> p = 0.2;  # error_probability

julia> recovery, success_prob = tn_decode(code, syndrome, p; contract_fn=mps_contract(8))
([0, 0, 0, 3, 3, 0, 1, 3, 0], 0.8250539267483997)

julia> find_syndrome(code, recovery) == syndrome
true
```
"""
function tn_decode(
    code::TensorNetworkCode, syndrome::AbstractVector{Int}, error_probability;
    contract_fn=basic_contract()
)
    pure_error = find_pure_error(code, syndrome)

    coset_probabilities = array(
        _contract(code, pure_error, error_probability, contract_fn)
    )

    if length(coset_probabilities) != 4
        error("decoder not set up for multiple logicals yet!")
    end

    # find best coset index
    coset_probabilities = abs.(coset_probabilities)
    best_coset_idx = argmax(coset_probabilities)

    # evaluate correction
    best_logical = best_coset_idx - 1  # logical in format 0=I, 1=X, 2=Y, 3=Z
    powers = Int[best_logical in (1, 2), best_logical in (2, 3)]  # convert to bsf
    correction = pauli_product_pow(code.logicals, powers)
    correction = pauli_product.(correction, pure_error)

    # evaluate success probability
    predicted_success_rate = coset_probabilities[best_coset_idx] / sum(coset_probabilities)

    return correction, predicted_success_rate
end

"""
    _contract(
        code, pure_error::AbstractVector{Int}, error_prob::Float64, contract_fn
    )

Contracts the tensor network described by `code_graph`.  Uses any specified `contract_fn`.
"""
function _contract(code, pure_error::AbstractVector{Int}, error_prob::Float64, contract_fn)
    L = Int64(sqrt(num_qubits(code)))
    imatrix = [ITensor(1) for i in 1:L, j in 1:L]

    for i in 1:L, j in 1:L
        node =  -(L*(i-1) + j)
        tensor = create_virtual_tensor(code, node)

        neighbours = physical_neighbours(code, node)
        for neighbour in neighbours
            indices = edge_indices(code, Set([node, neighbour]))
            α = findfirst(x->hastags(x, "physical"), indices)
            index = indices[α]
            tensor = tensor *
            physical_tensor(index, error_prob, pure_error[neighbour])
        end
        imatrix[i,j] = tensor
    end

    return contract_fn(imatrix)
end

end

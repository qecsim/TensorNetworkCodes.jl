module TNDecode

#imports
using ..TensorNetworkCodes: TensorNetworkCode
using ..TensorNetworkCodes: edges, edge_indices, node_indices, node_types
using ..TensorNetworkCodes: code_to_Itensor, create_virtual_tensor, physical_tensor
using ..TensorNetworkCodes: pauli_commutation, pauli_product, pauli_product_pow
using ..TensorNetworkCodes: num_qubits
using ..TensorNetworkCodes: find_pure_error, find_syndrome
using ..TensorNetworkCodes: random_pauli_error
using ITensors: ITensors  # only imported to avoid `contract` name clash
using ITensors: Index, ITensor, MPO, MPS, array, hastags
using Statistics: mean, stdm

#exports
export basic_contract, mps_contract, tn_decode
export TN_decoder
export compare_code_success_empirical, compare_code_success_predicted
export TN_decoding_success_prob

########## Should be in core/types.jl
function _physical_neighbours(code,node)

    all_edges = edges(code)
    neighbour_edges = filter(x->node in x,all_edges)
    neighbour_nodes = union(neighbour_edges...)

    return filter(x->x>0,neighbour_nodes)

end
##########

"""
    basic_contract

Given a matrix of ITensors (`imatrix`), contracts exactly working
across row by row.
"""
function basic_contract(imatrix; bond_dim = 10)
    output = ITensor(1)
    for i in 1:size(imatrix)[1], j in 1:size(imatrix)[2]
        output = output * imatrix[i,j]
    end

    return output
end

"""
    mps_contract

Given a matrix of ITensors (`imatrix`), contracts approximately by
treating the first row as an MPS and the subsequent rows as MPOs.
"""
function mps_contract(imatrix; bond_dim = 10)
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

"""
    _contract(
        code,
        pure_error::Array{Int64,1},
        error_prob::Float64;
        contraction_function=basic_contract,
        bond_dim = 10)

Contracts the tensor network described by `code_graph`.  Uses any
specified `contraction_function`.
"""
function _contract(
        code,
        pure_error::Array{Int64,1},
        error_prob::Float64;
        contraction_function=basic_contract,
        bond_dim = 10)

    L = Int64(sqrt(num_qubits(code)))
    imatrix = [ITensor(1) for i in 1:L, j in 1:L]

    for i in 1:L, j in 1:L
        node =  -(L*(i-1) + j)
        tensor = create_virtual_tensor(code,node)

        neighbours = _physical_neighbours(code,node)
        for neighbour in neighbours
            indices = edge_indices(code,Set([node,neighbour]))
            α = findfirst(x->hastags(x,"physical"),indices)
            index = indices[α]
            tensor = tensor *
            physical_tensor(index,error_prob,pure_error[neighbour])
        end
        imatrix[i,j] = tensor
    end

    return contraction_function(imatrix,bond_dim=bond_dim)

end

"""
    tn_decode(
        code::TensorNetworkCode, syndrome::AbstractVector{Int}, error_probability::Real;
        contract_fn=basic_contract, bond_dim=10
    ) -> Vector{Int}

TODO: functions that return contract function with parameters such as bond_dim.
TODO: doc
"""
function tn_decode(
    code::TensorNetworkCode, syndrome::AbstractVector{Int}, error_probability;
    contract_fn=basic_contract, bond_dim=10
)
    pure_error = find_pure_error(code,syndrome)

    coset_probabilities = array(_contract(
        code,
        pure_error,
        error_probability;
        contraction_function=contract_fn,
        bond_dim=bond_dim)
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
    TN_decoder(
        code,
        syndrome,
        error_prob;
        error_model = "depolarizing",
        contraction_function = simple_contract)

Decodes a `TensorNetworkCode` given the `syndrome` and a choice of
`TN_contraction_function`.
"""
function TN_decoder(
        code::TensorNetworkCode,
        syndrome::Array{Int64,1},
        error_prob::Float64;
        error_model = "depolarizing",
        contraction_function = basic_contract,
        bond_dim=10)

    # physical tensors
    pure_error = find_pure_error(code,syndrome)

    coset_probabilities = array(
        _contract(
        code,
        pure_error,
        error_prob;
        contraction_function=contraction_function,
        bond_dim = bond_dim)
        )

    if length(coset_probabilities) != 4
        error("decoder not set up for multiple logicals yet!")
    end

    coset_probabilities = abs.(coset_probabilities)

    best_coset = argmax(coset_probabilities) - 1
    powers = [0,0]
    if best_coset in [1,2]
        powers[1] = 1
    end
    if best_coset in [3,2]
        powers[2] = 1
    end

    correction = pauli_product_pow(code.logicals,powers)
    correction = pauli_product.(correction,pure_error)

    return correction
end

"""
    TN_decoding_success_prob(
        code,
        syndrome,
        error_prob;
        error_model = "depolarizing",
        contraction_function = basic_contract,
        bond_dim=10)

For a `TensorNetworkCode` given error `syndrome` and a choice of `contraction_function`,
finds the probability of successfully decoding.
"""
function TN_decoding_success_prob(
        code::TensorNetworkCode,
        syndrome::Array{Int64,1},
        error_prob::Float64;
        error_model = "depolarizing",
        contraction_function = basic_contract,
        bond_dim=10)

    # physical tensors
    pure_error = find_pure_error(code,syndrome)

    coset_probabilities = array(
        _contract(
        code,
        pure_error,
        error_prob;
        contraction_function=contraction_function,
        bond_dim = bond_dim)
        )

    coset_probabilities = abs.(coset_probabilities)
    best_coset = argmax(coset_probabilities)

    return coset_probabilities[best_coset]/(sum(coset_probabilities))
end

function compare_code_success_empirical(
        code1,
        code2,
        error_probabilities,
        num_iter;
        contraction_function = basic_contract,
        χ=10)


    P = length(error_probabilities)
    code1_success = zeros(P)
    code2_success = zeros(P)

    n = num_qubits(code1)
    if n != num_qubits(code2)
        error("different size codes")
    end

    for α in 1:P
        for _ in 1:num_iter
            error_prob = error_probabilities[α]

            error = random_pauli_error(n,error_prob)

            # code 1
            syndrome1 = find_syndrome(code1,error)

            code1_correction = TN_decoder(
                code1,
                syndrome1,
                error_prob;
                contraction_function = basic_contract,
                bond_dim=χ)

            effect_on_code = pauli_product.(error,code1_correction)
            if (pauli_commutation.(Ref(effect_on_code),code1.logicals)
                == zeros(Int64,length(code1.logicals)))
                code1_success[α] += 1
            end

            # code 2
            syndrome2 = find_syndrome(code2,error)

            code2_correction = TN_decoder(
                code2,
                syndrome2,
                error_prob;
                contraction_function = mps_contract,
                bond_dim=χ)

            effect_on_code = pauli_product.(error,code2_correction)
            if (pauli_commutation.(Ref(effect_on_code),code2.logicals)
                == zeros(Int64,length(code2.logicals)))
                code2_success[α] += 1
            end


        end
    end

    return [code1_success ./num_iter,code2_success ./num_iter]
end

function compare_code_success_predicted(
        code1,
        code2,
        error_probabilities,
        num_iter;
        contraction_function = basic_contract,
        χ=10)

    P = length(error_probabilities)
    code1_success = zeros(P)
    code1_stderr = zeros(P)
    code2_success = zeros(P)
    code2_stderr = zeros(P)

    n = num_qubits(code1)
    if n != num_qubits(code2)
        error("different size codes")
    end

    for α in 1:P
        data1 = zeros(num_iter)
        data2 = zeros(num_iter)

        for β in 1:num_iter
            error_prob = error_probabilities[α]
            error = random_pauli_error(n,error_prob)

            # code 1
            syndrome1 = find_syndrome(code1,error)

            code1_succ = TN_decoding_success_prob(
                code1,
                syndrome1,
                error_prob;
                error_model = "depolarizing",
                contraction_function = basic_contract,
                bond_dim=χ)

            data1[β] = code1_succ

            # code 2
            syndrome2 = find_syndrome(code2,error)

            code2_succ = TN_decoding_success_prob(
                code2,
                syndrome2,
                error_prob;
                error_model = "depolarizing",
                contraction_function = basic_contract,
                bond_dim=χ)

            data2[β] = code2_succ

        end

        code1_success[α] = mean(data1)
        code1_stderr[α] = stdm(data1,mean(data1))

        code2_success[α] = mean(data2)
        code2_stderr[α] = stdm(data2,mean(data2))

    end

    return [code1_success,code1_stderr ./ sqrt(num_iter),
        code2_success,code2_stderr ./ sqrt(num_iter)]
    # really it's only stderr (standard error) now after the division
end

end

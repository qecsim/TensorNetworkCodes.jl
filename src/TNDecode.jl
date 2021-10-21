module TNDecode

#imports
using ..TensorNetworkCodes: TNCode
using ..TensorNetworkCodes: edges, edge_indices, node_indices, node_types
using ..TensorNetworkCodes: code_to_Itensor
using ..TensorNetworkCodes: pauli_commutation, pauli_product, pauli_product_pow
using ..TensorNetworkCodes: num_qubits
using ..TensorNetworkCodes: get_syndrome, find_pure_error
using ..TensorNetworkCodes: random_pauli_error
using ITensors: ITensors  # only imported to avoid `contract` name clash
using ITensors: Index, ITensor, array, hastags
using Statistics: mean, stdm

#exports
export basic_contract, MPS_MPO_contract, TN_decoder
export compare_code_success_empirical, compare_code_success_predicted
export TN_decoding_success_prob

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
function create_virtual_tensor(code::TNCode,node::Int64)

    graph = code.code_graph

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

########## Should be in core/types.jl
function physical_neighbours(code,node)

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
    MPS_MPO_contract

Given a matrix of ITensors (`imatrix`), contracts approximately by
treating the first row as an MPS and the subsequent rows as MPOs.
"""
function MPS_MPO_contract(imatrix; bond_dim = 10)
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
    contract(
        code,
        pure_error::Array{Int64,1},
        error_prob::Float64;
        contraction_function=basic_contract,
        bond_dim = 10)

Contracts the tensor network described by `code_graph`.  Uses any
specified `contraction_function`.
"""
function contract(
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

        neighbours = physical_neighbours(code,node)
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
    TN_decoder(
        code,
        syndrome,
        error_prob;
        error_model = "depolarizing",
        contraction_function = simple_contract)

Decodes a `TNCode` given the `syndrome` and a choice of `TN_contraction_function`.
"""
function TN_decoder(
        code::TNCode,
        syndrome::Array{Int64,1},
        error_prob::Float64;
        error_model = "depolarizing",
        contraction_function = basic_contract,
        bond_dim=10)

    # physical tensors
    pure_error = find_pure_error(code,syndrome)

    coset_probabilities = array(
        contract(
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

For a `TNCode` given error `syndrome` and a choice of `contraction_function`,
finds the probability of successfully decoding.
"""
function TN_decoding_success_prob(
        code::TNCode,
        syndrome::Array{Int64,1},
        error_prob::Float64;
        error_model = "depolarizing",
        contraction_function = basic_contract,
        bond_dim=10)

    # physical tensors
    pure_error = find_pure_error(code,syndrome)

    coset_probabilities = array(
        contract(
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
            syndrome1 = get_syndrome(code1,error)

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
            syndrome2 = get_syndrome(code2,error)

            code2_correction = TN_decoder(
                code2,
                syndrome2,
                error_prob;
                contraction_function = MPS_MPO_contract,
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
            syndrome1 = get_syndrome(code1,error)

            code1_succ = TN_decoding_success_prob(
                code1,
                syndrome1,
                error_prob;
                error_model = "depolarizing",
                contraction_function = basic_contract,
                bond_dim=χ)

            data1[β] = code1_succ

            # code 2
            syndrome2 = get_syndrome(code2,error)

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

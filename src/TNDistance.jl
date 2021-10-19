module TNDistance

#imports
using ..TensorNetworkCodes: TNCode
using ..TensorNetworkCodes: all_cosets, code_to_Itensor, identity_coset
using ..TensorNetworkCodes: nodes, node_indices, node_types
using ..TensorNetworkCodes: pauli_weight
using ITensors: Index, ITensor, array, contract, dim, hastags, inds
using StatsPlots: groupedbar

#exports
export OperatorWeights
export operator_weights_plot
export TN_distance, TN_operator_weights

"""
    OperatorWeights

Stores operator weight distribution for stabilizers, all code operators
(stabilizers + all logical representatives), as well as distance.
"""
struct OperatorWeights
    stabilizer_weights::Array{Int64,1}
    all_operator_weights::Array{Int64,1}
    distance::Int64

    function OperatorWeights(stabilizer_weights,all_operator_weights)
        logical_weights = all_operator_weights - stabilizer_weights
        index = findfirst(x->x>0,logical_weights)

        distance =  index -1
        new(stabilizer_weights,all_operator_weights,distance)
    end
end
OperatorWeights() = OperatorWeights(Int64[0],Int64[1])

"""
    simple_weight_tensor(physical_indices,out_index)

Returns a tensor that, when contracted with a code tensor, gives the weights
of all Pauli operators described by that code tensor.
"""
function simple_weight_tensor(
        physical_indices::Array{Index{Int64},1},
        out_index::Index{Int64})

    all_indices = vcat(physical_indices,[out_index])
    dims = dim.(all_indices)
    array = fill(0,dims...)

    n = length(physical_indices)

    for γ in 0:4^n-1
        index = digits!(zeros(Int64,n),γ,base=4)
        w = pauli_weight(index)
        index = index .+ 1

        full_index = vcat(index,w+1)
        array[full_index...] = 1
    end

    return ITensor(array,all_indices...)
end

"""
    addition_tensor(in_index1,in_index2,out_index)

Gives a tensor that is one if `out_index` equals  the sum of `in_index1`
and `in_index2`.  Warning: since Julia is one indexed, zero is replaced by
1, so that the tensor is nonzero if `in_index1 + in_index2 - 1 = out_index`.
"""
function addition_tensor(in_index1,in_index2,out_index)

    in_dim1 = dim(in_index1)
    in_dim2 = dim(in_index2)
    out_dim = dim(out_index)
    dims = [in_dim1,in_dim2,out_dim]

    array = fill(0,dims...)

    for α in 1:in_dim1, β in 1:in_dim2, γ in 1:out_dim
        if α + β - 2 == γ - 1
            array[α,β,γ] = 1
        end
    end

    return ITensor(array,in_index1,in_index2,out_index)
end

function addition_tensor(in_index1,in_index2)

    in_dim1 = dim(in_index1)
    in_dim2 = dim(in_index2)
    out_dim = in_dim1 + in_dim2 -2 + 1  #(indices start from zero in Julia)
    out_index = Index(out_dim,"out_index")

    return addition_tensor(in_index1,in_index2,out_index)
end

"""
    TN_weights(code;cosets=all_cosets)

Returns a tensor describing the weights of all operators in the `cosets` of
the code.  (`cosets` is applied to each `seed_code` separately.)
"""
function TN_weights(
        code::TNCode;
        cosets=all_cosets,
        truncate_to=size(code)+1)

    code_tensors = ITensor[]

    for node in nodes(code)
        if node < 0
            code_type = node_types(code,node)
            seed_code = code.seed_codes[code_type]

            k = length(seed_code.logicals)/2
            logical_indices = [Index(4,"logical") for _ in 1:k]
            physical_indices = node_indices(code,node)
            code_tensor = code_to_Itensor(seed_code,logical_indices,physical_indices)

            code_tensor = cosets(code_tensor)

            push!(code_tensors,code_tensor)
        end
    end

    for α in 1:length(code_tensors)
        indices = inds(code_tensors[α])  #tuple
        physical_indices = filter(x->hastags(x,"physical"),indices)
        num_phys = length(physical_indices)
        output_index = Index(num_phys + 1,"weight")

        local_weight_tensor = simple_weight_tensor([physical_indices...],output_index)
        code_tensors[α] = code_tensors[α]*local_weight_tensor
    end


    indices = inds(code_tensors[1])
    β = findfirst(x->hastags(x,"weight"),indices)
    in_index1 = indices[β]

    for α in 2:length(code_tensors)
        indices = inds(code_tensors[α])
        β = findfirst(x->hastags(x,"weight"),indices)
        in_index2 = indices[β]

        in_dim1 = dim(in_index1)
        in_dim2 = dim(in_index2)
        total_out_dim = in_dim1 + in_dim2 -2 + 1

        if total_out_dim <= truncate_to
            out_dim = total_out_dim
        else
            out_dim = truncate_to
        end

        out_index = Index(out_dim,"outind")

        add_tensor = addition_tensor(in_index1,in_index2,out_index)
        code_tensors[α] = code_tensors[α] * add_tensor

        in_index1 = out_index
    end

    return contract(code_tensors)

end

"""
   TN_operator_weights(code)

Returns `OperatorWeights` which includes the number of stabilizers
of each weight, as well as the number of logical representatives
(excluding identity) of each weight, as well as the code distance.
"""
function TN_operator_weights(
        code::TNCode;
        truncate_to=size(code)+1)

    stabilizer_weights = array(TN_weights(
            code,
            cosets=identity_coset,
        truncate_to=truncate_to))

    all_operator_weights = array(TN_weights(
            code,
            cosets=all_cosets,
            truncate_to=truncate_to))

    return OperatorWeights(stabilizer_weights,all_operator_weights)
end

"""
   TN_distance(code)

Returns the code distance calculated by contracting a tensor network.
"""
function TN_distance(
        code::TNCode;
        truncate_to=size(code)+1)

    op_weights = TN_operator_weights(
        code::TNCode;
        truncate_to=truncate_to)

    return op_weights.distance
end

"""
    operator_weights_plot(weights)

Plots a bar plot of the number of operators of each weight.
This includes all code operators, all stabilizers, and all
non-identity logicals.
"""
function operator_weights_plot(
        weights::OperatorWeights;
        truncate_to = length(weights.stabilizer_weights))

    all_ops = weights.all_operator_weights[1:truncate_to]
    stabs = weights.stabilizer_weights[1:truncate_to]

    max_y_val = maximum(all_ops)

    logical_weights = all_ops - stabs

    # for log plot, so zero doesn't give -infinity
    logical_weights = logical_weights .+ 0.001
    all_ops = all_ops .+ 0.001
    stabs = stabs .+ 0.001

    bar_data = hcat(
        stabs,
        logical_weights)
#         all_ops)

    L = length(logical_weights)

    labels = repeat(
        ["Stabilizers", "Logicals"],
#             "Logicals + stabilizers"],
        inner = L)

    return groupedbar(
        repeat(collect(0:L-1),outer = 2),
        bar_data,
        bar_position = :dodge,
        group = labels,
        xlabel = "Weight",
        ylabel = "Number",
        yaxis=:log,
        title = "Number of operators of given weight",
        bar_width = 0.67,
        lw = 0,
        framestyle = :box,
        legend = :topleft,
        ylims = (0.1,10*max_y_val),
        xticks = (0:L-1))
end

end

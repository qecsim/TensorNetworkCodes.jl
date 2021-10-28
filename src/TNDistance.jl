module TNDistance

#imports
using ..TensorNetworkCodes: TensorNetworkCode
using ..TensorNetworkCodes: all_cosets, code_to_Itensor, identity_coset
using ..TensorNetworkCodes: nodes, node_indices, node_types
using ..TensorNetworkCodes: num_qubits
using ..TensorNetworkCodes: pauli_weight
using ITensors: Index, ITensor, array, contract, dim, hastags, inds
using StatsPlots: groupedbar

#exports
export OperatorWeights
export operator_weights_plot
export tn_distance, tn_operator_weights

"""
    OperatorWeights(stabilizer_weights,all_operator_weights)

Stores operator weight distribution for stabilizers and all code operators
(stabilizers + all logical representatives).  Also stores distance.  Fields:

    stabilizer_weights::Array{Int64,1}
    all_operator_weights::Array{Int64,1}
    distance::Int64

Note that since julia is 1-indexed `stabilizer_weights[j]` gives the number of
stabilizers of weight j-1.
"""
struct OperatorWeights
    stabilizer_weights::Array{Int64,1}
    all_operator_weights::Array{Int64,1}
    distance::Int64

    function OperatorWeights(stabilizer_weights,all_operator_weights)
        logical_weights = all_operator_weights - stabilizer_weights
        index = findfirst(x->x>0,logical_weights)

        if index == nothing
            distance = 0
        else
            distance =  index -1
        end
        new(stabilizer_weights,all_operator_weights,distance)
    end
end
OperatorWeights() = OperatorWeights(Int64[0],Int64[1])

"""
    simple_weight_tensor(physical_indices,out_index) -> ITensor

Returns a tensor that, when contracted with a code tensor, gives the weights
of all Pauli operators described by that code tensor.
"""
function _simple_weight_tensor(
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
    _addition_tensor(in_index1,in_index2,out_index)->ITensor

Gives a tensor that is one if `out_index` equals  the sum of `in_index1`
and `in_index2`.  Warning: since Julia is one indexed, zero is replaced by
1, so that the tensor is nonzero if `in_index1 + in_index2 - 1 = out_index`.
"""
function _addition_tensor(in_index1,in_index2,out_index)

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

function _addition_tensor(in_index1,in_index2)

    in_dim1 = dim(in_index1)
    in_dim2 = dim(in_index2)
    out_dim = in_dim1 + in_dim2 -2 + 1  #(indices start from zero in Julia)
    out_index = Index(out_dim,"out_index")

    return _addition_tensor(in_index1,in_index2,out_index)
end

"""
    _tn_weights(code::TensorNetworkCode;cosets=all_cosets,
    truncate_to=num_qubits(code)+1) -> ITensor

Returns a tensor describing the weights of all operators in the `cosets` (either
`all_cosets` or `identity_coset`) of
the code.  (`cosets` is applied to each `seed_code` separately.)
"""
function _tn_weights(
        code::TensorNetworkCode;
        cosets=all_cosets,
        truncate_to=num_qubits(code)+1)

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

        local_weight_tensor = _simple_weight_tensor([physical_indices...],output_index)
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

        add_tensor = _addition_tensor(in_index1,in_index2,out_index)
        code_tensors[α] = code_tensors[α] * add_tensor

        in_index1 = out_index
    end

    return contract(code_tensors)

end

"""
    tn_operator_weights(code::TensorNetworkCode;truncate_to=num_qubits(code)+1)
    -> OperatorWeights

Returns [`OperatorWeights`](@ref) which includes the number of stabilizers
of each weight, as well as the number of logical representatives
(excluding identity) of each weight, as well as the code distance.

# Examples
```jldoctest
julia> using TensorNetworkCodes.TNDistance

julia> code = TensorNetworkCode(five_qubit_code());

julia> code = contract(code,deepcopy(code),[[1,1],[3,3]]);

julia> tn_operator_weights(code)
OperatorWeights([1, 0, 0, 0, 9, 0, 6], [1, 0, 9, 24, 99, 72, 51], 2)
```
"""
function tn_operator_weights(
        code::TensorNetworkCode;
        truncate_to=num_qubits(code)+1)

    stabilizer_weights = array(_tn_weights(
            code,
            cosets=identity_coset,
        truncate_to=truncate_to))

    all_operator_weights = array(_tn_weights(
            code,
            cosets=all_cosets,
            truncate_to=truncate_to))

    return OperatorWeights(stabilizer_weights,all_operator_weights)
end

"""
    tn_distance(code::TensorNetworkCode;truncate_to=num_qubits(code)+1)
    -> Int64

Returns the code distance calculated by contracting a tensor network.

# Examples
```jldoctest
julia> using TensorNetworkCodes.TNDistance

julia> code = TensorNetworkCode(steane_code());

julia> code = contract(code,deepcopy(code),[[1,1],[3,3]]); # contract two copies

julia> tn_distance(code) # this code has poor distance!
2
```
"""
function tn_distance(
        code::TensorNetworkCode;
        truncate_to=num_qubits(code)+1)

    op_weights = tn_operator_weights(
        code::TensorNetworkCode;
        truncate_to=truncate_to)

    return op_weights.distance
end

"""
    operator_weights_plot(weights::OperatorWeights;truncate_to = length(weights.stabilizer_weights))

Plots a bar plot (on a log scale) of the number of operators of each weight.
This includes all code operators (logicals plus stabilizers) and all stabilizers plotted
separately.
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

    L = length(logical_weights)

    labels = repeat(
        ["Stabilizers", "Logicals"],
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

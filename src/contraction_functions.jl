# COMBINE FUNCTIONS

"""
    combine(code1::SimpleCode, code2::SimpleCode) -> SimpleCode
    combine(code1::TensorNetworkCode, code2::TensorNetworkCode) -> TensorNetworkCode

Return a new code that is the tensor product of the given codes. Physically equivalent to
preparing two codes on different sets of physical qubits.

# Examples
```jldoctest
julia> code = combine(five_qubit_code(), steane_code());  # combine 5 and 7 qubit codes

julia> num_qubits(code)
12
```
"""
function combine(code1::SimpleCode, code2::SimpleCode)
    n1 = num_qubits(code1)
    n2 = num_qubits(code2)

    stabilizers = Vector{Int}[]
    logicals = Vector{Int}[]
    pure_errors = Vector{Int}[]

    for stabilizer in code1.stabilizers
        push!(stabilizers, vcat(stabilizer, zeros(Int, n2)))
    end
    for stabilizer in code2.stabilizers
        push!(stabilizers, vcat(zeros(Int, n1), stabilizer))
    end

    for logical in code1.logicals
        push!(logicals, vcat(logical, zeros(Int, n2)))
    end
    for logical in code2.logicals
        push!(logicals, vcat(zeros(Int, n1), logical))
    end

    for pure_error in code1.pure_errors
        push!(pure_errors, vcat(pure_error, zeros(Int, n2)))
    end
    for pure_error in code2.pure_errors
        push!(pure_errors, vcat(zeros(Int, n1), pure_error))
    end

    return SimpleCode("", stabilizers, logicals, pure_errors)
end
function combine(code1::TensorNetworkCode, code2::TensorNetworkCode)
    # code2 = (code1 === code2) ? deepcopy(code2) : code2  # copy if codes identical
    
    # copy & change indices of one of the codes.  necessary to avoid
    # duplicate ITensor indices
    if num_qubits(code1) < num_qubits(code2)
        code1 = new_indices(code1)
    else 
        code2 = new_indices(code2)
    end

    new_code = combine(SimpleCode(code1), SimpleCode(code2))
    new_code_graph = _combine(code1.code_graph, code2.code_graph)
    new_seed_codes = merge(code1.seed_codes, code2.seed_codes)

    return TensorNetworkCode(new_code, new_code_graph, new_seed_codes)
end

"""
    _combine(graph1,graph2)

Combines two `CodeGraphs` with physical and virtual nodes and
preserves their metadata (`type`, `indices`, `coords` and 'qubit'
labels).
"""
function _combine(code_graph1::CodeGraph, code_graph2::CodeGraph)
    n1 = num_nodes(code_graph1)
    Labels = nodes(code_graph1)
    n_qubits1 = length(filter(x -> x > 0, Labels))

    new_code_graph2 = _shift_keys(code_graph2, n_qubits1, n1)

    coords = merge(code_graph1.coords, new_code_graph2.coords)
    node_types = merge(code_graph1.node_types, new_code_graph2.node_types)
    edge_types = merge(code_graph1.edge_types, new_code_graph2.edge_types)
    node_indices = merge(code_graph1.node_indices, new_code_graph2.node_indices)
    edge_indices = merge(code_graph1.edge_indices, new_code_graph2.edge_indices)

    return CodeGraph(coords, node_types, edge_types, node_indices, edge_indices)
end

"""
    _shift_keys(code_graph,n_qubits,n_vert)

When merging two `Graphs`, the second `Graph` has all its edges
relabelled.  So to merge two `CodeGraphs`, it is necessary to preemptively
shift the keys of all the dictionaries in the second `CodeGraph` to
allow for this.

`code_graph` is the second of the two `CodeGraphs` to be merged.  `n_vert`
and `n` are the number of vertices and qubits respectively of the first
`CodeGraph` to be merged.
"""
function _shift_keys(code_graph::CodeGraph, n_qubits::Int64, n_vert::Int64)

    new_coords = Dict{Int64,Vector{Float64}}()
    new_node_types = Dict{Int64,String}()
    new_edge_types = Dict{Set{Int64},String}()
    new_node_indices = Dict{Int64,Vector{Index{Int64}}}()
    new_edge_indices = Dict{Set{Int64},Vector{Index{Int64}}}()

    n_virtual = n_vert - n_qubits

    # Update node data
    for (key, value) in coords(code_graph)
        if key > 0
            new_coords[key + n_qubits] = value
        else
            new_coords[key - n_virtual] = value
        end
    end
    for (key, value) in node_types(code_graph)
        if key > 0
            new_node_types[key + n_qubits] = value
        else
            new_node_types[key - n_virtual] = value
        end
    end
    for (key, value) in node_indices(code_graph)
        if key < 0
            new_node_indices[key - n_virtual] = value
        end
    end

    # Edge data is a bit trickier
    for old_edge in edges(code_graph)
        new_edge = old_edge .+ 0 # converts to vector

        for α in 1:2
            if new_edge[α] > 0
                new_edge[α] = new_edge[α] + n_qubits
            else
                new_edge[α] = new_edge[α] - n_virtual
            end
        end

        new_edge = Set(new_edge)

        new_edge_types[new_edge] = edge_types(code_graph, old_edge)
        new_edge_indices[new_edge] = edge_indices(code_graph, old_edge)
    end

    return CodeGraph(
        new_coords,
        new_node_types,
        new_edge_types,
        new_node_indices,
        new_edge_indices
    )
end

# FUSION FUNCTIONS

"""
    fusion(code::SimpleCode, qubit_pair::AbstractVector{Int}) -> SimpleCode
    fusion(code::SimpleCode, qubit_pairs) -> SimpleCode

    fusion(code::TensorNetworkCode, qubit_pair::AbstractVector{Int}) -> TensorNetworkCode
    fusion(code::TensorNetworkCode, qubit_pairs) -> TensorNetworkCode

Return a new code that results from fusing the physical qubits with labels given by each
qubit pair. The new code has two fewer physical qubits, for each qubit pair, but the same
number of logical qubits. Physically equivalent to updating stabilizers, logicals and pure
errors after measuring ``XX`` and ``ZZ`` on each pair of qubits.

The versions that take `qubit_pairs` take iterables of `AbstractVector{Int}`. An
`ErrorException` is thrown if the fusion is not possible, i.e., if the logical
degrees of freedom would not be preserved by the measurements.

# Examples
```jldoctest
julia> code = steane_code();

julia> num_qubits(code), length(code.logicals) ÷ 2
(7, 1)

julia> fused_code = fusion(code, [1, 2]);

julia> num_qubits(fused_code), length(fused_code.logicals) ÷ 2
(5, 1)
```
"""
function fusion(code::SimpleCode, qubit_pair::AbstractVector{Int})
    stabilizers = deepcopy(code.stabilizers)
    pure_errors = deepcopy(code.pure_errors)
    logicals = deepcopy(code.logicals)

    annoying_ops = _annoying_stabilizers(stabilizers, qubit_pair)
    if length(annoying_ops) != 0
        @debug "There were annoying operators!"
    end
    for logical in logicals
        if !pauli_are_commuting(vcat(annoying_ops, [logical]))
            error("Fusion is not possible between selected qubit pair!")
        end
    end

    useful_inds = _find_useful_stabilizer_indices(stabilizers, qubit_pair)
    if length(useful_inds) != 0
        sort!(useful_inds)
        useful_ops = [stabilizers[ind] for ind in useful_inds]
        deleteat!(stabilizers, sort!(useful_inds))
        deleteat!(pure_errors, sort!(useful_inds))

        stabilizers = _make_ready(stabilizers, useful_ops, qubit_pair)
        logicals = _make_ready(logicals, useful_ops, qubit_pair)
        pure_errors = _make_ready(pure_errors, useful_ops, qubit_pair)
    end

    _remove_qubits!(stabilizers, qubit_pair)
    _remove_qubits!(pure_errors, qubit_pair)
    if length(logicals) != 0
        _remove_qubits!(logicals, qubit_pair)
    end

    if length(useful_inds) != 2
        new_stabilizers = Array{Int64,1}[]
        new_pure_errors = Array{Int64,1}[]
        for α in 1:length(stabilizers)
            if pauli_are_independent(vcat(new_stabilizers, [stabilizers[α]]))
                push!(new_stabilizers, stabilizers[α])
                push!(new_pure_errors, pure_errors[α])
            end
        end
        pure_errors = new_pure_errors
        stabilizers = new_stabilizers
    end

    return SimpleCode(" ", stabilizers, logicals, pure_errors)
end
function fusion(code::SimpleCode, qubit_pairs)
    output_code = deepcopy(code)
    qubit_pair_list = collect(qubit_pairs) # convert interable to list
    # qubit labels change after each fusion, so must account for this:
    _update_qubit_pairs!(qubit_pair_list)

    for qubit_pair in qubit_pair_list
        if num_qubits(output_code)  == 0
            return SimpleCode()
        end
        output_code = fusion(output_code, qubit_pair)
    end

    return output_code
end
function fusion(code::TensorNetworkCode, qubit_pairs)
    new_code = fusion(SimpleCode(code), qubit_pairs)
    new_code_graph = _fusion(code.code_graph, qubit_pairs)
    return TensorNetworkCode(new_code, new_code_graph, code.seed_codes)
end
function fusion(code::TensorNetworkCode, qubit_pair::AbstractVector{Int})
    return fusion(code, [qubit_pair])
end

"""
    _fusion(code_graph,qubit_pair)

Returns a `CodeGraph` after fusion is performed on the qubits in
`qubit_pair`.
"""
function _fusion(code_graph::CodeGraph, qubit_pair::AbstractVector{Int})
    new_edge_types = deepcopy(code_graph.edge_types)
    new_node_indices = deepcopy(code_graph.node_indices)
    new_edge_indices = deepcopy(code_graph.edge_indices)

    # Find new edge
    edges_to_remove = edges(code_graph)
    edges_to_remove =
    filter(x -> length(intersect(qubit_pair, x)) != 0, edges_to_remove)

    virtual1 = setdiff(edges_to_remove[1], qubit_pair)
    virtual2 = setdiff(edges_to_remove[2], qubit_pair)
    new_edge = union(virtual1, virtual2)

    if length(new_edge) == 1
        println("Self contraction occurred.  Tensor-network decoding & distance algorithms may not work!")
        @goto skip_new_edge_stuff
    end

    # Add edge data
    new_edge_types[new_edge] = "bond"
    new_index = Index(4, "bond")
    new_edge_indices[new_edge] = [new_index]

    # replace node indices with new bond index
    indices_to_replace = [new_edge_indices[edge][1] for edge in edges_to_remove]
    for label in new_edge # this is cool: new_edge is a Set
        index_list = new_node_indices[label]
        for α in 1:length(index_list)
            if index_list[α] in indices_to_replace
                index_list[α] = new_index
            end
        end
    end

    @label skip_new_edge_stuff

    output = CodeGraph(
        code_graph.coords,
        code_graph.node_types,
        new_edge_types,
        new_node_indices,
        new_edge_indices
    )

    sort!(qubit_pair)
    qubit_pair[2] = qubit_pair[2] - 1
    for qubit in qubit_pair
        output = _remove_node(output, qubit)
    end

    return output
end
function _fusion(code_graph::CodeGraph, qubit_pairs::Array{Array{Int64,1},1})
    output_graph = code_graph
    for qubit_pair in qubit_pairs
        output_graph = _fusion(output_graph, qubit_pair)
    end

    return output_graph
end

function _annoying_stabilizers(
    stabilizers::AbstractVector{<:AbstractVector{Int}},
    qubit_pair::AbstractVector{Int}
)
    i, j = qubit_pair
    n = length(stabilizers[1])
    output = Array{Array{Int64,1},1}()

    for σ in [1,2,3]

        for α in 1:length(stabilizers)
            if !pauli_are_commuting([[σ,σ],[stabilizers[α][i],stabilizers[α][j]]])
                @goto not_this_one
            end
        end

        operator = zeros(Int64, n)
        operator[i] = σ
        operator[j] = σ
        push!(output, operator)

        if length(output) == 2
            break
        end

        @label not_this_one
    end

    return output
end

"""
    _find_useful_stabilizer_indices(stabilizers,qubit_pair)

Finds two `stabilizers` that can be used to make other operators
ready for fusion.  Basically it returns two operators that are independent
from each other and XX, YY, ZZ on the qubits in `qubit_pair`.
"""
function _find_useful_stabilizer_indices(
    stabilizers::AbstractVector{<:AbstractVector{Int}},
    qubit_pair::AbstractVector{Int}
)
    indices = Int64[]

    for α in 1:length(stabilizers)
        operator = stabilizers[α]
        if operator[qubit_pair[1]] == operator[qubit_pair[2]]
            continue
        end

        if length(indices) == 0
            push!(indices, α)
            continue
        end

        if length(indices) == 1
            useful_operator = stabilizers[indices[1]]
            new_operator = pauli_product.(operator, useful_operator)

            if new_operator[qubit_pair[1]] != new_operator[qubit_pair[2]]
                push!(indices, α)
                return indices
            end
        end
    end

    return indices
end

"""
    _make_ready(operator,useful_ops,qubit_pair)

Makes an `operator` ready for fusion on `qubit_pair`, i.e.,
multiplies a combination of operators in `useful_ops` until the
resulting operator has the same Paulis on both qubits in `qubit_pair`.
"""
function _make_ready(
    operator::AbstractVector{Int},
    useful_ops::AbstractVector{<:AbstractVector{Int}},
    qubit_pair::AbstractVector{Int}
)
    if _is_ready(operator, qubit_pair)
        return operator
    end

    for useful_op in useful_ops
        new_operator = pauli_product.(operator, useful_op)
        if _is_ready(new_operator, qubit_pair)
            return new_operator
        end
    end

    if length(useful_ops) == 2
        operator2 = pauli_product.(useful_ops[1], useful_ops[2])
        new_operator = pauli_product.(operator, operator2)
        if _is_ready(new_operator, qubit_pair)
            return new_operator
        end
    end
end
_make_ready(
    operators::AbstractVector{<:AbstractVector{Int}},
    useful_ops::Array{Array{Int64,1}},
    qubit_pair::AbstractVector{Int}
) = _make_ready.(operators, Ref(useful_ops), Ref(qubit_pair))

"""
    _is_ready(operator,qubit_pair)

Checks if `operator` is ready for fusion on `qubit_pair`, i.e.,
does `operator` have the same Paulis on both qubits.
"""
function _is_ready(operator::AbstractVector{Int}, qubit_pair::AbstractVector{Int})
    return (operator[qubit_pair[1]] == operator[qubit_pair[2]])
end

"""
    _remove_qubits!(operators, qubits_to_remove)

Given a set of operators, remove entries corresponding to qubits being removed.
"""
function _remove_qubits!(
    operators::AbstractVector{<:AbstractVector{Int}},
    qubits_to_remove::AbstractVector{Int}
)
    for α in 1:length(operators)
        deleteat!(operators[α], sort!(qubits_to_remove))
    end
end

function _update_qubit_pairs!(qubit_pairs::AbstractVector{<:AbstractVector{Int}})
    num_pairs = length(qubit_pairs)
    for α in 1:num_pairs
        for x in 1:2
            for β in 1:α - 1
                if qubit_pairs[β][1] <= qubit_pairs[α][x]
                    qubit_pairs[α][x] -= 1
                end
                if qubit_pairs[β][2] <= qubit_pairs[α][x]
                    qubit_pairs[α][x] -= 1
                end
            end
        end
    end
end

"""
    _remove_node(code_graph,label)

Removes a node from a `CodeGraph`; returns a new `CodeGraph`.
"""
function _remove_node(code_graph::CodeGraph, label::Int)

    # Remove graph vertex and update values in `nodes`
#     # A little annoying because of how rem_vertex! works
#     # in LightGraphs.jl
#     g = code_graph.graph
#     end_label = nodes_to_labels(code_graph)[end]
#     removed_node = nodes(code_graph,label)

#     new_g = deepcopy(g)
#     rem_vertex!(new_g,removed_node)

    # Update dictionaries with nodes as keys
    old_labels = nodes(code_graph)
    new_labels = Dict{Int64,Int64}()
    for old_label in old_labels
        if old_label < label
            new_labels[old_label] = old_label
        elseif old_label > label
           new_labels[old_label] = old_label - 1
        end
    end

    new_coords = _update_keys(code_graph.coords, new_labels)
    new_node_types = _update_keys(code_graph.node_types, new_labels)
    filter!(x -> x.first < 0, new_labels)
    new_node_indices = _update_keys(code_graph.node_indices, new_labels)

    # Update dictionaries with edges as keys
    old_edges = edges(code_graph)
    new_edges = Dict{Set{Int64},Set{Int64}}()
    for old_edge in old_edges

        new_edge = old_edge .+ 0 # converts to vector

        if label in new_edge
            continue
        end
        for α in 1:2
            if new_edge[α] > label
                new_edge[α] = new_edge[α] - 1
            end
        end

        new_edge = Set(new_edge)
        new_edges[old_edge] = new_edge
    end

    new_edge_types = _update_keys(code_graph.edge_types, new_edges)
    new_edge_indices = _update_keys(code_graph.edge_indices, new_edges)

    return CodeGraph(
        new_coords,
        new_node_types,
        new_edge_types,
        new_node_indices,
        new_edge_indices
    )
end

# maybe don't bother including types
function _update_keys(
    dict::Dict{Int,T}, new_keys::Dict{Int,Int}
) where T <: Union{Int,Vector{Float64},String,Vector{Index{Int}}}

    value_type = typeof(dict[findfirst(x -> true, dict)])  # this is lame
    new_dict = Dict{Int64,value_type}()

    for old_key in keys(new_keys)
        new_key = new_keys[old_key]
        new_dict[new_key] = dict[old_key]
    end

    return new_dict
end
function _update_keys(
    dict::Dict{Set{Int},T}, new_keys::Dict{Set{Int},Set{Int}}
) where T <: Union{Int,Vector{Float64},String,Vector{Index{Int}}}

    value_type = typeof(dict[findfirst(x -> true, dict)])  # this is lame
    new_dict = Dict{Set{Int64},value_type}()

    for old_key in keys(new_keys)
        new_key = new_keys[old_key]
        new_dict[new_key] = dict[old_key]
    end

    return new_dict
end

# CONTRACT FUNCTIONS

"""
    contract(
        code1::TensorNetworkCode,
        code2::TensorNetworkCode,
        qubit_pair::AbstractVector{Int}
    ) -> TensorNetworkCode

    contract(
        code1::TensorNetworkCode,
        code2::TensorNetworkCode,
        qubit_pairs
    ) -> TensorNetworkCode

Return a new code that results from combining the codes and fusing physical qubits
identified by the qubit pairs. The first and second elements of a qubit pair refer to qubit
labels from the first and second codes, respectively.

This is equivalent to [`combine`](@ref) followed by [`fusion`](@ref), with the qubit pair
labels referring to the code qubit labels before combining. The version that takes
`qubit_pairs` takes iterables of `AbstractVector{Int}`. An `ErrorException` is thrown if the
fusion is not possible.

See also: [`combine`](@ref), [`fusion`](@ref).

# Examples
```jldoctest
julia> code1 = TensorNetworkCode(five_qubit_code());

julia> code2 = TensorNetworkCode(steane_code());

julia> contracted_code = contract(code1, code2, [[1, 2], [2, 7]]);

julia> num_qubits(contracted_code), length(contracted_code.logicals) ÷ 2
(8, 2)

julia> fusion(combine(code1, code2), [[1, 7], [2, 12]]);  # equivalent code
```
"""
function contract(code1::TensorNetworkCode, code2::TensorNetworkCode, qubit_pairs)
    code2 = (code1 === code2) ? deepcopy(code2) : code2  # copy if codes identical

    n = num_qubits(code1)
    qubit_pair_list = [[pair[1], pair[2] + n] for pair in qubit_pairs]

    output_code = combine(code1, code2)
    output_code = fusion(output_code, qubit_pair_list)
    return output_code
end
function contract(
    code1::TensorNetworkCode,
    code2::TensorNetworkCode,
    qubit_pair::AbstractVector{Int}
)
    return contract(code1, code2, [qubit_pair])
end

"""
    contract_by_coords(code1::TensorNetworkCode,code2::TensorNetworkCode)
        -> TensorNetworkCode

Return a new code that results from combining the codes and fusing physical qubits with
coincident coordinates.

This is equivalent to [`combine`](@ref) followed by [`fusion`](@ref), with the fusion qubit
pairs being those with coincident coordinates. An `ErrorException` is thrown if the fusion
is not possible.

See also: [`combine`](@ref), [`fusion`](@ref).

# Examples
```jldoctest
julia> code1 = TensorNetworkCode(five_qubit_code());

julia> code2 = TensorNetworkCode(steane_code());

julia> set_coords!(code1, 1, [5, 5]); set_coords!(code2, 2, [5, 5]); # qubits 1 & 2 coincide

julia> set_coords!(code1, 2, [6, 6]); set_coords!(code2, 7, [6, 6]); # qubits 2 & 7 coincide

julia> contracted_code = contract_by_coords(code1, code2);

julia> num_qubits(contracted_code), length(contracted_code.logicals) ÷ 2
(8, 2)

julia> contract(code1, code2, [[1, 2], [2, 7]]);  # equivalent code
```
"""
function contract_by_coords(code1::TensorNetworkCode, code2::TensorNetworkCode)
    graph1 = code1.code_graph
    graph2 = code2.code_graph
    qubit_pairs = Vector{Int}[]

    for label1 in nodes(graph1), label2 in nodes(graph2)
        if (graph1.node_types[label1] != "physical" ||
            graph2.node_types[label2] != "physical")
            continue
        end

        coord1 = coords(graph1, label1)
        coord2 = coords(graph2, label2)
        if coord1 ≈ coord2
            push!(qubit_pairs, [label1,label2])
        end
    end

    if length(qubit_pairs) == 0
        return combine(code1, code2)
    end

    return contract(code1, code2, qubit_pairs)
end

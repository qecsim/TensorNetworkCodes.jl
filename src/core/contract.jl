"""
    remove_qubits!(operators,qubits_to_remove)

Given a set of operators, remove entries corresponding to qubits
being removed.
"""
function remove_qubits!(operators::Array{Array{Int64,1}},
        qubits_to_remove::Array{Int64,1})

    for α in 1:length(operators)
        deleteat!(operators[α],sort!(qubits_to_remove))
    end
end





"""
    merge(code1,code2)

Merges two `SimpleCodes`.  Physically equivalent to preparing two
codes on different sets of qubits.  Mathematically, this is the
tensor product of two codes.
"""
function Base.merge(code1::SimpleCode,code2::SimpleCode)

    n1 = size(code1)
    n2 = size(code2)

    stabilizers = Array{Int64,1}[]
    logicals = Array{Int64,1}[]
    pure_errors = Array{Int64,1}[]

    for stabilizer in code1.stabilizers
        push!(stabilizers,vcat(stabilizer,zeros(Int,n2)))
    end
    for stabilizer in code2.stabilizers
        push!(stabilizers,vcat(zeros(Int,n1),stabilizer))
    end

    for logical in code1.logicals
        push!(logicals,vcat(logical,zeros(Int,n2)))
    end
    for logical in code2.logicals
        push!(logicals,vcat(zeros(Int,n1),logical))
    end

    for pure_error in code1.pure_errors
        push!(pure_errors,vcat(pure_error,zeros(Int,n2)))
    end
    for pure_error in code2.pure_errors
        push!(pure_errors,vcat(zeros(Int,n1),pure_error))
    end


    return SimpleCode("",stabilizers,logicals,pure_errors)
end





"""
    is_ready(operator,qubit_pair)

Checks if `operator` is ready for fusion on `qubit_pair`, i.e.,
does `operator` have the same Paulis on both qubits.
"""
function is_ready(operator::Array{Int64,1},qubit_pair::Array{Int64,1})
    if operator[qubit_pair[1]] == operator[qubit_pair[2]]
        return true
    else
        return false
    end
end





"""
    make_ready(operator,useful_ops,qubit_pair)

Makes an `operator` ready for fusion on `qubit_pair`, i.e.,
multiplies a combination of operators in `useful_ops` until the
resulting operator has the same Paulis on both qubits in `qubit_pair`.
"""
function make_ready(
        operator::Array{Int64,1},
        useful_ops::Array{Array{Int64,1}},
        qubit_pair::Array{Int64,1})

    if is_ready(operator,qubit_pair)
        return operator
    end

    for useful_op in useful_ops
        new_operator = pauli_product.(operator,useful_op)
        if is_ready(new_operator,qubit_pair)
            return new_operator
        end
    end

    if length(useful_ops) == 2
        operator2 = pauli_product.(useful_ops[1],useful_ops[2])
        new_operator = pauli_product.(operator,operator2)
        if is_ready(new_operator,qubit_pair)
            return new_operator
        end
    end
end



make_ready(operators::Array{Array{Int64,1},1},useful_ops::Array{Array{Int64,1}},
        qubit_pair::Array{Int64,1}) =
make_ready.(operators,Ref(useful_ops),Ref(qubit_pair))





"""
    find_useful_stabilizer_indices(stabilizers,qubit_pair)

Finds two `stabilizers` that can be used to make other operators
ready for fusion.  Basically it returns two operators that are independent
from each other and XX, YY, ZZ on the qubits in `qubit_pair`.
"""
function find_useful_stabilizer_indices(stabilizers::Array{Array{Int64,1},1},
        qubit_pair::Array{Int64,1})

    indices = Int64[]

    for α in 1:length(stabilizers)
        operator = stabilizers[α]
        if operator[qubit_pair[1]] == operator[qubit_pair[2]]
            continue
        end

        if length(indices) == 0
            push!(indices,α)
            continue
        end

        if length(indices) == 1
            useful_operator = stabilizers[indices[1]]
            new_operator = pauli_product.(operator,useful_operator)

            if new_operator[qubit_pair[1]] != new_operator[qubit_pair[2]]
                push!(indices,α)
                return indices
            end
        end
    end

    return indices
end





function update_qubit_pairs!(code1::QuantumCode,qubit_pairs::Array{Array{Int64,1},1})

    num_pairs = length(qubit_pairs)
    for α in 1:num_pairs
        for x in 1:2
            for β in 1:α-1
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





function annoying_stabilizers(
        stabilizers::Array{Array{Int64,1},1},
        qubit_pair::Array{Int64,1})

    i,j = qubit_pair
    n = length(stabilizers[1])
    output = Array{Array{Int64,1},1}()

    for σ in [1,2,3]
        operator = zeros(Int64,n)
        operator[i] = σ
        operator[j] = σ

        if pauli_are_commuting(vcat(stabilizers,[operator]))
            push!(output,operator)
        end
        if length(output) == 2
            break
        end
    end

    return output
end





"""
    fusion(code,qubit_pair)

Physically equivalent to updating stabilizers, logicals and pure_errors
after measuring XX and ZZ on the qubits in `qubit_pair`.

Returns a `SimpleCode` with two fewer physical qubits but the same number
of logical qubits.
"""
function fusion(code::SimpleCode,qubit_pair::Array{Int64,1})

    stabilizers = deepcopy(code.stabilizers)
    pure_errors = deepcopy(code.pure_errors)
    logicals = deepcopy(code.logicals)


    annoying_ops = annoying_stabilizers(stabilizers,qubit_pair)
    if length(annoying_ops) != 0
        println("There were annoying operators!")
    end
    for logical in logicals
        if !pauli_are_commuting(vcat(annoying_ops,[logical]))
            println("Can't contract!")
            return SimpleCode()
        end
    end


    useful_inds = find_useful_stabilizer_indices(stabilizers,qubit_pair)
    if length(useful_inds) != 0
        sort!(useful_inds)
        useful_ops = [stabilizers[ind] for ind in useful_inds]
        deleteat!(stabilizers,sort!(useful_inds))
        deleteat!(pure_errors,sort!(useful_inds))

        stabilizers = make_ready(stabilizers,useful_ops,qubit_pair)
        logicals = make_ready(logicals,useful_ops,qubit_pair)
#         pure_errors = make_ready(pure_errors,useful_ops,qubit_pair)
#         fix_pure_errors!(pure_errors,stabilizers)
    end


    remove_qubits!(stabilizers,qubit_pair)
    if length(logicals) != 0
        remove_qubits!(logicals,qubit_pair)
    end


    if length(useful_inds) != 2
        new_stabilizers = Array{Int64,1}[]
        for α in 1:length(stabilizers)
            if are_they_independent(vcat(new_stabilizers,[stabilizers[α]]))
                push!(new_stabilizers,stabilizers[α])
            end
        end
        stabilizers = new_stabilizers
    end


#     remove_qubits!(pure_errors,qubit_pair)
    pure_errors = generate_pure_errors(stabilizers)

    return SimpleCode(" ",stabilizers,logicals,pure_errors)
end



function fusion(code::SimpleCode,qubit_pairs::Array{Array{Int64,1},1})

    output_code = deepcopy(code)
    # qubit labels change after each fusion, so must account
    # for this:
    update_qubit_pairs!(code,qubit_pairs)

    for qubit_pair in qubit_pairs
        if size(output_code)  == 0
            return SimpleCode()
        end
        output_code = fusion(output_code,qubit_pair)
    end

    return output_code
end









# maybe don't bother including types
function update_keys(
        dict::Dict{Int64,T},
        new_keys::Dict{Int64,Int64}) where T <: Union{
    Int64,
    Vector{Float64},
    String,Vector{Index{Int64}}}

    value_type = typeof(dict[findfirst(x->true,dict)])  # this is lame
    new_dict = Dict{Int64,value_type}()

    for old_key in keys(new_keys)
        new_key = new_keys[old_key]
        new_dict[new_key] = dict[old_key]
    end

    return new_dict
end


function update_keys(
        dict::Dict{Set{Int64},T},
        new_keys::Dict{Set{Int64},Set{Int64}}) where T <: Union{
    Int64,
    Vector{Float64},
    String,Vector{Index{Int64}}}

    value_type = typeof(dict[findfirst(x->true,dict)])  # this is lame
    new_dict = Dict{Set{Int64},value_type}()

    for old_key in keys(new_keys)
        new_key = new_keys[old_key]
        new_dict[new_key] = dict[old_key]
    end

    return new_dict
end





"""
    shift_keys(code_graph,n_qubits,n_vert)

When merging two `Graphs`, the second `Graph` has all its edges
relabelled.  So to merge two `CodeGraphs`, it is necessary to preemptively
shift the keys of all the dictionaries in the second `CodeGraph` to
allow for this.

`code_graph` is the second of the two `CodeGraphs` to be merged.  `n_vert`
and `n` are the number of vertices and qubits respectively of the first
`CodeGraph` to be merged.
"""
function shift_keys(code_graph::CodeGraph,n_qubits::Int64,n_vert::Int64)

    new_coords = Dict{Int64,Vector{Float64}}()
    new_node_types = Dict{Int64,String}()
    new_edge_types = Dict{Set{Int64},String}()
    new_node_indices = Dict{Int64,Vector{Index{Int64}} }()
    new_edge_indices = Dict{Set{Int64},Vector{Index{Int64}} }()


    n_virtual = n_vert - n_qubits


    # Update node data
    for (key,value) in coords(code_graph)
        if key > 0
            new_coords[key + n_qubits] = value
        else
            new_coords[key - n_virtual] = value
        end
    end
    for (key,value) in node_types(code_graph)
        if key > 0
            new_node_types[key + n_qubits] = value
        else
            new_node_types[key - n_virtual] = value
        end
    end
    for (key,value) in node_indices(code_graph)
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

        new_edge_types[new_edge] = edge_types(code_graph,old_edge)
        new_edge_indices[new_edge] = edge_indices(code_graph,old_edge)
    end


    return CodeGraph(
        new_coords,
        new_node_types,
        new_edge_types,
        new_node_indices,
        new_edge_indices)
end





"""
    merge(graph1,graph2)

Merges two `CodeGraphs` with physical and virtual nodes and
preserves their metadata (`type`, `indices`, `coords` and 'qubit'
labels).
"""
function Base.merge(code_graph1::CodeGraph,code_graph2::CodeGraph)

    n1 = num_nodes(code_graph1)
    Labels = nodes(code_graph1)
    n_qubits1 = length(filter(x->x>0,Labels))

    new_code_graph2 = shift_keys(code_graph2,n_qubits1,n1)

    coords = merge(code_graph1.coords,new_code_graph2.coords)
    node_types = merge(code_graph1.node_types,new_code_graph2.node_types)
    edge_types = merge(code_graph1.edge_types,new_code_graph2.edge_types)
    node_indices = merge(code_graph1.node_indices,new_code_graph2.node_indices)
    edge_indices = merge(code_graph1.edge_indices,new_code_graph2.edge_indices)


    return CodeGraph(
        coords,
        node_types,
        edge_types,
        node_indices,
        edge_indices)
end





"""
    rem_node(code_graph,label)

Removes a node from a `CodeGraph`; returns a new `CodeGraph`.
"""
function rem_node(code_graph::CodeGraph,label::Int64)

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


    new_coords = update_keys(code_graph.coords,new_labels)
    new_node_types = update_keys(code_graph.node_types,new_labels)
    filter!(x->x.first<0,new_labels)
    new_node_indices = update_keys(code_graph.node_indices,new_labels)


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

    new_edge_types = update_keys(code_graph.edge_types,new_edges)
    new_edge_indices = update_keys(code_graph.edge_indices,new_edges)


    return CodeGraph(
    new_coords,
    new_node_types,
    new_edge_types,
    new_node_indices,
    new_edge_indices)
end





"""
    fusion(code_graph,qubit_pair)

Returns a `CodeGraph` after fusion is performed on the qubits in
`qubit_pair`.
"""
function fusion(code_graph::CodeGraph,qubit_pair::Array{Int64,1})

    new_edge_types = deepcopy(code_graph.edge_types)
    new_node_indices = deepcopy(code_graph.node_indices)
    new_edge_indices = deepcopy(code_graph.edge_indices)


    # Find new edge
    edges_to_remove = edges(code_graph)
    edges_to_remove =
    filter(x->length(intersect(qubit_pair,x))!=0,edges_to_remove)

    virtual1 = setdiff(edges_to_remove[1],qubit_pair)
    virtual2 = setdiff(edges_to_remove[2],qubit_pair)
    new_edge = union(virtual1,virtual2)

    if length(new_edge) == 1
        println("Self contraction occurred.  Contraction algorithms may not work!")
        @goto skip_new_edge_stuff
    end

    # Add edge data
    new_edge_types[new_edge] = "bond"
    new_index = Index(4,"bond")
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
    new_edge_indices)


    sort!(qubit_pair)
    qubit_pair[2] = qubit_pair[2]-1
    for qubit in qubit_pair
        output = rem_node(output,qubit)
    end


    return output
end



# Another method with multiple `qubit_pairs`.
function fusion(code_graph::CodeGraph,qubit_pairs::Array{Array{Int64,1},1})

    output_graph = code_graph
    for qubit_pair in qubit_pairs
        output_graph = fusion(output_graph,qubit_pair)
    end

    return output_graph
end








"""
    merge(code1,code2) -> typeof(code1)

Takes the (disjoint) union of two error correcting codes.
"""
function Base.merge(code1::TNCode,code2::TNCode)

    new_code = merge(SimpleCode(code1),SimpleCode(code2))
    new_code_graph = merge(code1.code_graph,code2.code_graph)
    new_seed_codes = merge(code1.seed_codes,code2.seed_codes)

    return TNCode(new_code,new_code_graph,new_seed_codes)
end





"""
    fusion(code,qubit_pairs)

If possible, implements fusion on pairs of qubits in `qubit_pairs` to return
a new `QuantumCode` of the same type as the input.
"""
function fusion(code::TNCode,qubit_pairs::Array{Array{Int64,1},1})

    new_code = fusion(SimpleCode(code),qubit_pairs)
    new_code_graph = fusion(code.code_graph,qubit_pairs)

    return TNCode(new_code,new_code_graph,code.seed_codes)
end



fusion(code::TNCode,qubit_pair::Array{Int64,1}) =
fusion(code,[qubit_pair])





"""
    contract(code1,code2,vertex_pairs) -> TNCode

Combines two `TN_codes` to give a new one by contracting each pair of indices
corresponding to qubits in `qubit_pairs`.
"""
function contract(
        code1::TNCode,
        code2::TNCode,
        qubit_pairs::Array{Array{Int64,1},1})

    qubit_pairs = [[qubit_pairs[m][1],qubit_pairs[m][2]+size(code1)]
            for m in 1:length(qubit_pairs)]

    output_code = merge(code1,code2)
    output_code = fusion(output_code,qubit_pairs)

end





"""
    combine_by_coordinates(code1::TNCode,code2::TNCode)

Combines two `TNCode` by finding out which physical qubits live on
overlapping coordinates and contracting those.
"""
function combine_by_coordinates(code1::TNCode,code2::TNCode)

    graph1 = code1.code_graph
    graph2 = code2.code_graph
    qubit_pairs = Array{Int64,1}[]

    for label1 in nodes(graph1), label2 in nodes(graph2)
        if (graph1.node_types[label1] != "physical" ||
            graph2.node_types[label2] != "physical")
            continue
        end

        coord1 = coords(graph1,label1)
        coord2 = coords(graph2,label2)
        if coord1 ≈ coord2
            push!(qubit_pairs,[label1,label2])
        end
    end


    if length(qubit_pairs) == 0
       return merge(code1,code2)
    end

    return contract(code1,code2,qubit_pairs)
end

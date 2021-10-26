

"""
    num_nodes(code::TensorNetworkCode)

Returns how many nodes are in the `TensorNetworkCode`, which is the sum
of the number of qubits and number of virtual tensors.
"""
num_nodes(code_graph::CodeGraph) = length(code_graph.coords)



num_nodes(code::TensorNetworkCode) = num_nodes(code.code_graph)





# Useful lookup functions
for data in [:coords,:node_types, :node_indices,
        :edge_types, :edge_indices]
    eval(quote
#         access(code::TensorNetworkCode,data::Symbol) = code.code_graph.$data
        $data(code_graph::CodeGraph) = code_graph.$data
        $data(code::TensorNetworkCode) = code.code_graph.$data
    end)
end



for data = (:coords, :node_types, :node_indices)
    eval(quote
        $data(code_graph::CodeGraph,label::Int64) = code_graph.$data[label]
        $data(code::TensorNetworkCode,label::Int64) = code.code_graph.$data[label]
    end)
end



for data = (:edge_types, :edge_indices)
    eval(quote
        $data(code_graph::CodeGraph,label::Set{Int64}) = code_graph.$data[label]
        $data(code::TensorNetworkCode,label::Set{Int64}) = code.code_graph.$data[label]
    end)
end





# Useful modification functions
for data in [:coords,:node_types, :node_indices,
        :edge_types, :edge_indices]
    function_name = Symbol(:set_,data,:!)
    eval(quote
        function $function_name(code::TensorNetworkCode,new_data)
            code.code_graph.$data = new_data
            return nothing
        end
        function $function_name(code_graph::CodeGraph,new_data)
            code_graph.$data = new_data
            return nothing
        end
    end)
end



for data = (:coords, :node_types, :node_indices)
    function_name = Symbol(:set_,data,:!)
    eval(quote
        function $function_name(code::TensorNetworkCode,label::Int64,new_data)
            code.code_graph.$data[label] = new_data
            return nothing
        end
        function $function_name(code_graph::CodeGraph,label::Int64,new_data)
            code_graph.$data[label] = new_data
            return nothing
        end
    end)
end



for data = (:edge_types, :edge_indices)
    function_name = Symbol(:set_,data,:!)
    eval(quote
        function $function_name(code::TensorNetworkCode,label::Set{Int64},new_data)
            code.code_graph.$data[label] = new_data
            return nothing
        end
        function $function_name(code_graph::CodeGraph,label::Set{Int64},new_data)
            code_graph.$data[label] = new_data
            return nothing
        end
    end)
end





"""
    nodes(code)

Returns an ordered list of nodes corresponding to qubits and
virtual nodes for a `TensorNetworkCode`.
"""
function nodes(code_graph::CodeGraph)
    output = collect(keys(coords(code_graph)))
    sort!(output)
    return output
end



nodes(code::TensorNetworkCode) = nodes(code.code_graph)





"""
    edges(code)

Returns a list of `CodeEdges` for a `TensorNetworkCode`.
"""
function edges(code_graph::CodeGraph)
    output = collect(keys(code_graph.edge_types))
    return output
end



edges(code::TensorNetworkCode) = edges(code.code_graph)
















"""
    change_code_coords!(code, coords)

Assigns coordinates to each node of a `TN_code` including both physical
and virtual nodes.
"""
function set_coords!(
        code::TensorNetworkCode,
        new_coords::Vector{Vector{T}}) where T <: Real

    if num_nodes(code) != length(new_coords)
        error("number of nodes and coordinates don't match!")
    end

    Ω = nodes(code)
    for α in 1:length(Ω)
        ω = Ω[α]
        set_coords!(code,ω,new_coords[α])
    end
end





"""
    shift_coords!(code::TN_code, shift::Array{Int64,1})

Shifts coordinates of all nodes of a `TN_code`.
"""
function shift_coords!(
        code::TensorNetworkCode,
        shift::Array{T,1}) where T <: Real

    for node in nodes(code)
        new_coords = coords(code,node) + shift
        set_coords!(code,node,new_coords)
    end
end

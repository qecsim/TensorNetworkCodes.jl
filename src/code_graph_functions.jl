# getters, e.g. coords(code::TensorNetworkCode)
for data in [:coords, :node_types, :node_indices, :edge_types, :edge_indices]
    eval(quote
        @doc """
            $($data)(code::TensorNetworkCode)

        Return `$($data)` field of the code's [`CodeGraph`](@ref).
        """
        $data(code::TensorNetworkCode) = code.code_graph.$data
        $data(code_graph::CodeGraph) = code_graph.$data
    end)
end

# getters, e.g. coords(code::TensorNetworkCode, label::Int)
for data = (:coords, :node_types, :node_indices)
    eval(quote
        @doc """
            $($data)(code::TensorNetworkCode, label::Int)

        Return `$($data)[label]` for the node from the code's [`CodeGraph`](@ref).
        """
        $data(code::TensorNetworkCode, label::Int) = code.code_graph.$data[label]
        $data(code_graph::CodeGraph, label::Int) = code_graph.$data[label]
    end)
end

# getters, e.g. edge_types(code::TensorNetworkCode, label::Set{Int})
for data = (:edge_types, :edge_indices)
    eval(quote
        @doc """
            $($data)(code::TensorNetworkCode, label::Set{Int})

        Return `$($data)[label]` for the edge from the code's [`CodeGraph`](@ref).
        """
        $data(code::TensorNetworkCode, label::Set{Int}) = code.code_graph.$data[label]
        $data(code_graph::CodeGraph, label::Set{Int}) = code_graph.$data[label]
    end)
end

# NOTE: THESE SETTERS CANNOT WORK BECAUSE CodeGraph IS IMMUTABLE
# # setters, e.g. set_coords!(code::TensorNetworkCode, new_data)
# for data in [:coords, :node_types, :node_indices, :edge_types, :edge_indices]
#     function_name = Symbol(:set_,data,:!)
#     eval(quote
#         function $function_name(code::TensorNetworkCode, new_data)
#             code.code_graph.$data = new_data
#             return nothing
#         end
#         function $function_name(code_graph::CodeGraph, new_data)
#             code_graph.$data = new_data
#             return nothing
#         end
#     end)
# end

# setters, e.g. set_coords!(code::TensorNetworkCode, label:Int, new_data)
for data = (:coords, :node_types, :node_indices)
    function_name = Symbol(:set_,data,:!)
    eval(quote
        @doc """
            $($function_name)(code::TensorNetworkCode, label::Int, new_data)

        Set `$($data)[label]` for the node from the code's [`CodeGraph`](@ref).
        """
        function $function_name(code::TensorNetworkCode, label::Int, new_data)
            code.code_graph.$data[label] = new_data
            return nothing
        end
        function $function_name(code_graph::CodeGraph, label::Int, new_data)
            code_graph.$data[label] = new_data
            return nothing
        end
    end)
end

# setters, e.g. set_edge_types!(code::TensorNetworkCode, label:Set{Int}, new_data)
for data = (:edge_types, :edge_indices)
    function_name = Symbol(:set_,data,:!)
    eval(quote
        @doc """
            $($function_name)(code::TensorNetworkCode, label::Set{Int}, new_data)

        Set `$($data)[label]` for the edge from the code's [`CodeGraph`](@ref).
        """
        function $function_name(code::TensorNetworkCode, label::Set{Int}, new_data)
            code.code_graph.$data[label] = new_data
            return nothing
        end
        function $function_name(code_graph::CodeGraph, label::Set{Int}, new_data)
            code_graph.$data[label] = new_data
            return nothing
        end
    end)
end

function edges(code_graph::CodeGraph)
    output = collect(keys(code_graph.edge_types))
    return output
end
"""
    edges(code) -> Vector{Set{Int}}

Returns a list of edge labels for the code, corresponding to edges between nodes.
"""
edges(code::TensorNetworkCode) = edges(code.code_graph)

function nodes(code_graph::CodeGraph)
    output = collect(keys(coords(code_graph)))
    sort!(output)
    return output
end
"""
    nodes(code::TensorNetworkCode) -> Vector{Int}

Returns an ordered list of node labels for the code, corresponding to the qubits and virtual
nodes.
"""
nodes(code::TensorNetworkCode) = nodes(code.code_graph)

num_nodes(code_graph::CodeGraph) = length(code_graph.coords)
"""
    num_nodes(code::TensorNetworkCode) -> Int

Returns the number of nodes in the code, which is the sum of the number of qubits and the
number of virtual nodes.
"""
num_nodes(code::TensorNetworkCode) = num_nodes(code.code_graph)

"""
    physical_neighbours(code::TensorNetworkCode, node::Int) -> Set{Int}

Returns the node labels of physical qubits connected by an edge to the given node.
"""
function physical_neighbours(code,node)
    all_edges = edges(code)
    neighbour_edges = filter(x -> node in x, all_edges)
    neighbour_nodes = union(neighbour_edges...)
    return filter(x -> x > 0 && x != node, neighbour_nodes)
end

"""
    set_coords!(
        code::TensorNetworkCode,
        new_coords::AbstractVector{<:AbstractVector{<:Real}}
    )

Assign new coordinates to all (physical and virtual) nodes of the code, where the order of
the coordinates matches that of the nodes returned by [`nodes`](@ref), i.e. in ascending
order.
"""
function set_coords!(
        code::TensorNetworkCode,
        new_coords::AbstractVector{<:AbstractVector{<:Real}}
)
    num_nodes(code) == length(new_coords) || error("different number of nodes and coords!")
    for (i, node) in enumerate(nodes(code))
        set_coords!(code, node, new_coords[i])
    end
end

"""
    shift_coords!(code::TensorNetworkCode, shift::AbstractVector{<:Real})

Shift the coordinates of all (physical and virtual) nodes of the code by the given shift.
"""
function shift_coords!(code::TensorNetworkCode, shift::AbstractVector{<:Real})
    for node in nodes(code)
        new_coords = coords(code,node) + shift
        set_coords!(code, node, new_coords)
    end
end


"""
    new_indices(code::TensorNetworkCode) -> TensorNetworkCode

Deepcopies the code but changes all `ITensor` indices.
"""
function new_indices(code::TensorNetworkCode)
    new_code = deepcopy(code)

    for edge in edges(new_code)
        old_indices = edge_indices(new_code,edge)
        # Create new indices with same dims and tags
        new_indices = [Index(dim(old_index),tags(old_index)) for old_index in old_indices]
        for node in nodes(new_code)
            if node < 0
                new_node_indices = node_indices(new_code,node)
                for α in 1:length(old_indices)
                    β = findfirst(x->x==old_indices[α],new_node_indices)
                    if β == nothing
                        continue
                    end
                    new_node_indices[β] = new_indices[α]
                end
                set_node_indices!(new_code,node,new_node_indices)
            end
        end
        set_edge_indices!(new_code,edge,new_indices)
    end

    return new_code
end

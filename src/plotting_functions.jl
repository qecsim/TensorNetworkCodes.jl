"""
    plot_code(code::TensorNetworkCode; use_coords=true)

Return a plot, as returned by `GraphPlot.gplot`, of the code. If used in a Jupyter
notebook the plot is rendered after the cell.

Physical qubit nodes are coloured red and virtual tensors are coloured green.

# Examples
```julia
julia> code = TensorNetworkCode(five_qubit_code());

julia> plot_code(code)
```
![code plot](./plotting_functions/plot_code.svg)
"""
function plot_code(code::TensorNetworkCode; use_coords=true)
    Nodes = nodes(code.code_graph)
    graph = _light_graph(code)

    # nodelabels = [(node_types(code,node) == "physical") ? node : " "
    #     for node in Nodes]

    # nodecolours = [(node_types(code,node) == "physical") ?
    #     "tomato" : "lightgreen"
    #     for node in Nodes]

    nodelabels = fill("",length(Nodes))
    nodecolours = fill("",length(Nodes))

    for α in 1:length(Nodes)
        node = Nodes[α]
        if node_types(code,node) == "physical"
            nodelabels[α] = string(node)
            nodecolours[α] = "tomato"
        else
            code_name = node_types(code,node)
            seed_code = code.seed_codes[code_name]
            if length(seed_code.logicals) == 0
                nodecolours[α] = "lightgreen"
            else
                nodecolours[α] = "dodgerblue"
            end
        end
    end

    locs_x = Float64[]
    locs_y = Float64[]
    for node in Nodes
        locs = coords(code,node)
        push!(locs_x,locs[1])
        push!(locs_y,locs[2])
    end

    if use_coords
        return gplot(graph,locs_x,locs_y,nodefillc=nodecolours,nodelabel=nodelabels)
        # display(gplot(graph,locs_x,locs_y,nodefillc=nodecolours,nodelabel=nodelabels))
    else
        return gplot(graph,nodefillc=nodecolours,nodelabel=nodelabels)
        # display(gplot(graph,nodefillc=nodecolours,nodelabel=nodelabels))
    end
end

"""
    plot_operator(code::TensorNetworkCode, operator::AbstractVector{Int}; use_coords=true)

Return a plot, as returned by `GraphPlot.gplot`, of the operator. If used in a Jupyter
notebook the plot is rendered after the cell.

# Examples
```julia
julia> code = TensorNetworkCode(steane_code());

julia> plot_operator(code, code.stabilizers[1])
```
![operator plot](./plotting_functions/plot_operator.svg)
"""
function plot_operator(
    code::TensorNetworkCode, operator::AbstractVector{Int};
    use_coords = true
)
    Nodes = nodes(code.code_graph)
    graph = _light_graph(code)

    nodelabels = []
    nodecolours = String[]
    for node in Nodes
        node_type = node_types(code,node)

        if node_type == "physical"
            pauli = pauli_rep_change(operator[node])
            push!(nodelabels,pauli)
            if pauli == 'I'
                push!(nodecolours,"goldenrod1")
            else
                push!(nodecolours,"tomato")
            end
        else
            seed_code = code.seed_codes[node_type]
            if length(seed_code.logicals) == 0
                push!(nodecolours,"lightgreen")
            else
                push!(nodecolours,"dodgerblue")
            end
            push!(nodelabels," ")
        end
    end

    locs_x = Float64[]
    locs_y = Float64[]
    for node in Nodes
        locs = coords(code,node)
        push!(locs_x,locs[1])
        push!(locs_y,locs[2])
    end

    if use_coords
        return gplot(graph,locs_x,locs_y,nodefillc=nodecolours,nodelabel=nodelabels)
        # display(gplot(graph,locs_x,locs_y,nodefillc=nodecolours,nodelabel=nodelabels))
    else
        return gplot(graph,nodefillc=nodecolours,nodelabel=nodelabels)
        # display(gplot(graph,nodefillc=nodecolours,nodelabel=nodelabels))
    end
end

"""
    _light_graph(code_graph::CodeGraph)

Given a `CodeGraph`, returns a `LightGraph` for plotting.
"""
function _light_graph(code_graph::CodeGraph)

    Nodes = nodes(code_graph)
    virtual_nodes = filter(x->x<0,Nodes)
    n_virtual = length(virtual_nodes)

    Edges = edges(code_graph)

    output = Graph(length(Nodes))
    for edge in Edges
        graph_edge = _edge_to_graph_edge(edge,n_virtual)
        add_edge!(output,graph_edge...)
    end

    return output
end
_light_graph(code::TensorNetworkCode) = _light_graph(code.code_graph)

"""
    _edge_to_graph_edge(edge::Set{Int64},n_virtual::Int64)

Because `LightGraph` indices start at 1, but `CodeGraph` indices
include negative numbers, this converts `CodeGraph` edges to have
positive values.
"""
function _edge_to_graph_edge(edge::Set{Int64},n_virtual::Int64)

    output = [edge...] .+ n_virtual
    for α in 1:2
        if output[α] < n_virtual
            output[α] += 1
        end
    end

    return output
end

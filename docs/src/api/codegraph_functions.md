# Code graph functions

Functions to access and manipulate the [`CodeGraph`](@ref) associated with a
[`TensorNetworkCode`](@ref).

## Introductory notes

* Node labels are of type `Int` with negative and positive integers representing
    virtual nodes and physical qubits, respectively.
* Edge labels are of type `Set{Int}` and contain two elements corresponding to the
    labels of the nodes that the edge links.
* Node coordinates are of type `Vector{<:Real}` and contain two elements corresponding
    to 2-dimensional Cartesian coordinates.
* Node types are of type `String` and can take the values: "physical", ... TODO
* Edge types are of type `String` and can take the values: "physical", "bond", ... TODO
* Node and edge indices are of type `ITensor.Index`, ... TODO

## Specialized functions
```@docs
edges
nodes
num_nodes
set_coords!(::TensorNetworkCode, ::AbstractVector{<:AbstractVector{<:Real}})
shift_coords!(::TensorNetworkCode, ::AbstractVector{<:Real})
```

## Generic getters / setters
```@docs
coords
set_coords!(::TensorNetworkCode, ::Int, ::Any)
edge_indices
set_edge_indices!
edge_types
set_edge_types!
node_indices
set_node_indices!
node_types
set_node_types!
```

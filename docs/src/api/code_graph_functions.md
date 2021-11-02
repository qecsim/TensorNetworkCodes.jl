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
* Node types are of type `String` and can take the values: "physical" for physical qubits or else they take the name of a seed code for virtual nodes, e.g., "Five qubit code".
* Edge types are of type `String` and can take the values: "physical" for an edge between a virtual node and physical (i.e., qubit) node or "bond" for an edge between two virtual nodes, e.g., as occurs after contracting a "Five qubit code" and a "Steane code".
* Edge indices are of type `Vector{ITensors.Index{Int64}}}`.  These are lists of the indices necessary for tensor-network contraction (e.g., for decoding).
* Node indices are also of type `Vector{ITensors.Index{Int64}}}`.  These indices match the indices on the incident edges to that node.  This is necessary for, e.g., tensor-network decoding: when a tensor is created at a node, the first index in the list of node indices corresponds to the first qubit of the seed code, and so on.

## Specialized functions
```@docs
edges
nodes
num_nodes
physical_neighbours
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
new_indices
```

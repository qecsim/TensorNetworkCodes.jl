# Code graph functions

Functions to access and manipulate the [`CodeGraph`](@ref) associated with a
[`TensorNetworkCode`](@ref).

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

# Types

Types representing quantum error correcting codes. An abstract type, [`QuantumCode`](@ref),
defines the contract for all codes. [`SimpleCode`](@ref) is a basic implementation that can
be used as a seed code in a [`TensorNetworkCode`](@ref), which has an associated
[`CodeGraph`](@ref).

## Quantum code

```@docs
QuantumCode
```

## Simple code

```@docs
SimpleCode
SimpleCode()
SimpleCode(name, stabilizers, logicals)
SimpleCode(::TensorNetworkCode)
```

## Tensor-network code

```@docs
CodeGraph
TensorNetworkCode
TensorNetworkCode(::SimpleCode)
```

# Code functions

[`QuantumCode`](@ref) functions for verifying, evaluating, transforming and contracting.

## Basic

Functions for code properties or verification.

```@docs
num_qubits
verify_code
```

## Evaluation

Functions to evaluate operators or syndromes.

```@docs
find_distance_logicals
find_pure_error
find_pure_errors
find_syndrome
```

## Transformation

Functions to gauge, permute and purify codes.

```@docs
gauge
permute
purify
```

## Contraction

Functions to contract codes using the primitives of combining codes and fusing physical
qubits.

```@docs
combine
contract
contract_by_coords
fusion
```
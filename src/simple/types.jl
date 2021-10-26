"""
Abstract supertype for quantum error correcting codes.

Subtypes should include the following fields:

    stabilizers::Vector{Vector{Int}}
    logicals::Vector{Vector{Int}}
    pure_errors::Vector{Vector{Int}}

where `stabilizers` are a linearly independent collection of stabilizer generators,
`logicals` are representatives of logical operators ordered as ``X``-type, ``Z``-type, ...
for each logical qubit, and `pure_errors` (or destabilizers) are operators the anticommute
with the corresponding stabilizer generator and commute with all other stabilizer
generators. Each of `stabilizers`, `logicals` and `pure_errors` are collections of
multi-qubit Pauli operators, see [Pauli functions](@ref). The relationships between the
operators can be checked using [`verify_code`](@ref).

!!! note

    Transformation and contraction functions do not modify `stabilizers`, `logicals` or
    `pure_errors` fields directly but rather return new codes. Therefore these fields can be
    treated as immutable.

"""
abstract type QuantumCode end

"""
    SimpleCode <: QuantumCode

    SimpleCode(name, stabilizers, logicals, pure_errors)

Basic named implementation of a quantum code.

Fields in addition to [`QuantumCode`](@ref):

    name::String

This code contains no information on physical qubit layout but it can be used as a seed code
in a [`TensorNetworkCode`](@ref).
"""
struct SimpleCode <: QuantumCode
    name::String
    stabilizers::Vector{Vector{Int}}
    logicals::Vector{Vector{Int}}
    pure_errors::Vector{Vector{Int}}
end

"""
    SimpleCode()

Construct empty simple code with empty string name.
"""
SimpleCode() = SimpleCode("", Vector{Int}[], Vector{Int}[], Vector{Int}[])

"""
    SimpleCode(name, stabilizers, logicals)

Construct a simple code with pure errors evaluated using [`find_pure_errors`](@ref).

An `ErrorException` is thrown if the pure errors cannot be evaluated; for example, if the
stabilizers are not linearly independent.

```jldoctest
julia> code = SimpleCode(
           "5-qubit code",
           [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]],
           [[1, 1, 1, 1, 1], [3, 3, 3, 3, 3]]
       );

julia> code.pure_errors
4-element Vector{Vector{Int64}}:
 [0, 1, 0, 0, 0]
 [1, 3, 0, 0, 0]
 [3, 1, 0, 0, 0]
 [1, 0, 0, 0, 0]

julia> verify_code(code)
true
```
"""
function SimpleCode(name, stabilizers, logicals)
    # convert stabilizers for find_pure_errors
    typed_stabilizers = convert(Vector{Vector{Int}}, stabilizers)
    return SimpleCode(name, stabilizers, logicals, find_pure_errors(typed_stabilizers))
end

# REFACTORED ABOVE

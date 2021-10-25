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
"""
abstract type QuantumCode end

# REFACTORED ABOVE

"""
    `SimpleCode`

Type for simple quantum error correcting codes,
where we don't care about physical qubit layout.
"""
struct SimpleCode <: QuantumCode
    name::String
    stabilizers::Array{Array{Int64,1},1}
    logicals::Array{Array{Int64,1},1}
    pure_errors::Array{Array{Int64,1},1}
end



SimpleCode() = SimpleCode(
    "",
    Array{Int64,1}[],
    Array{Int64,1}[],
    Array{Int64,1}[])

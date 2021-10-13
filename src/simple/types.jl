"""
    `QuantumCode`

Abstract type for quantum error correcting codes.
"""
abstract type QuantumCode
#     should have at least these attributes:
#     stabilizers::Array{Array{Int64,1},1}
#     logicals::Array{Array{Int64,1},1}
#     pure_errors::Array{Array{Int64,1},1}
end





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

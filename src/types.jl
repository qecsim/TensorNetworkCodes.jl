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

Simple named implementation of a quantum code.

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

Construct an empty simple code with empty string name.
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

"""
    CodeGraph

Type aggregating all geometrical data used by [`TensorNetworkCode`](@ref), including
coordinates and `ITensor` indices.

See also [`TensorNetworkCode`](@ref).

!!! note

    `CodeGraph` is not typically accessed directly but rather created, read and updated
    using [Transformation](@ref), [Contraction](@ref) and [Code graph functions](@ref).
"""
struct CodeGraph
    coords::Dict{Int,Vector{T}} where T <: Real
    node_types::Dict{Int,String}
    edge_types::Dict{Set{Int},String}
    node_indices::Dict{Int,Vector{Index{Int}}}
    edge_indices::Dict{Set{Int},Vector{Index{Int}}}
end

"""
    TensorNetworkCode <: QuantumCode

    TensorNetworkCode(stabilizers, logicals, pure_errors, code_graph, seed_codes)

Tensor-network implementation of a quantum code, with associated code graph and seed codes.

Fields in addition to [`QuantumCode`](@ref):

    code_graph::CodeGraph
    seed_codes::Dict{String,SimpleCode}

The `code_graph` contains information on physical qubit layout. The `seed_codes` are simple
codes that are used to build the tensor-network code.

See also [`CodeGraph`](@ref).

!!! note

    The `code_graph` and `seed_codes` fields are not typically accessed directly but rather
    read and updated using [Transformation](@ref), [Contraction](@ref) and
    [Code graph functions](@ref).
"""
struct TensorNetworkCode <: QuantumCode
    stabilizers::Vector{Vector{Int}}
    logicals::Vector{Vector{Int}}
    pure_errors::Vector{Vector{Int}}
    code_graph::CodeGraph
    seed_codes::Dict{String,SimpleCode}
end

# non-public constructor
function TensorNetworkCode(
    code::SimpleCode,
    code_graph::CodeGraph,
    seed_codes::Dict{String,SimpleCode}
)
    return TensorNetworkCode(
        code.stabilizers, code.logicals, code.pure_errors, code_graph, seed_codes
    )
end

"""
    TensorNetworkCode(code::SimpleCode)

Construct a tensor-network code from a seed simple code and initialize the code graph with
with appropriate defaults.

The constructed code shares the `stabilizers`, `logicals` and `pure_errors` with the given
code.

See also [`SimpleCode`](@ref).

```jldoctest
julia> simple_code = five_qubit_code();

julia> tn_code = TensorNetworkCode(simple_code);

julia> verify_code(tn_code)
true
```
"""
function TensorNetworkCode(code::SimpleCode)
    n = num_qubits(code)
    # initialize code graph dictionaries
    coords = Dict{Int,Vector{Float64}}()
    node_types = Dict{Int,String}()
    edge_types = Dict{Set{Int},String}()
    node_indices = Dict{Int,Vector{Index{Int}}}()
    edge_indices = Dict{Set{Int},Vector{Index{Int}}}()
    # assign coords to qubits
    coords[-1] = [0.0, 0.0]
    vec = [0.0, -1.0]
    θ = 2 * pi / n
    R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
    for node in 1:n
        coords[node] = vec
        vec = R * vec
    end
    # label nodes
    node_types[-1] = code.name
    for label in 1:n
        node_types[label] = "physical"
    end
    # label edges
    for node in 1:n
        edge = Set([-1,node])
        edge_types[edge] = "physical"
    end
    # assign indices
    indices = [Index(4, "physical") for _ in 1:n]
    node_indices[-1] = indices
    for label in 1:n
        edge = Set([-1,label])
        edge_indices[edge] = [indices[label]]
    end
    # create code_graph
    code_graph = CodeGraph(coords, node_types, edge_types, node_indices, edge_indices)
    # return new code
    return TensorNetworkCode(
        code.stabilizers,
        code.logicals,
        code.pure_errors,
        code_graph,
        Dict([(code.name, code)])
    )
end

"""
     SimpleCode(code::TensorNetworkCode)

Construct a simple code from a tensor-network code with an empty string name.

The constructed code shares the `stabilizers`, `logicals` and `pure_errors` with the given
code.

See also [`TensorNetworkCode`](@ref).
"""
function SimpleCode(code::TensorNetworkCode)
    #TODO: add name="" as parameter?
    return SimpleCode("", code.stabilizers, code.logicals, code.pure_errors)
end

"""
Adaptors for using tensor-network codes within the
[Qecsim](https://github.com/dkt29/Qecsim.jl) framework.

Code, [`QecsimTNCode`](@ref), and decoder, [`QecsimTNDecoder`](@ref), implementations are
provided to be used with Qecsim error models for QEC simulations using the Qecsim `App`
module.
"""
module QecsimAdaptors

#imports
using ..TensorNetworkCodes: TensorNetworkCode, num_qubits
using ..TNDecode: basic_contract, mps_contract, tn_decode
using Qecsim

#exports
export QecsimTNCode, QecsimTNDecoder

"""
    QecsimTNCode <: Qecsim.StabilizerCode

    QecsimTNCode(code::TensorNetworkCode; distance=missing, label=nothing)

Qecsim stabilizer code implementation based on a tensor-network code.

Pauli operators such as `stabilizers` and `logicals` are converted to Qecsim's format on
construction. These are available via Qecsim API methods of the same name. The Qecsim `nkd`
method returns calculated values for `n` (number of physical qubits) and `k` (number of
logical qubits) and the given `distance` for `d`, if provided, otherwise `missing`. The
Qecsim `label` method returns the given `label`, if provided, otherwise the label is
constructed from `[n,k,d]`.

Public fields:

    tn_code::TensorNetworkCode  # the given tensor-network code

# Examples
```jldoctest
julia> using TensorNetworkCodes.QecsimAdaptors

julia> using Qecsim: label, nkd, validate

julia> tn_code = TensorNetworkCode(five_qubit_code());

julia> qs_code = QecsimTNCode(tn_code; distance=3);

julia> label(qs_code)
"QecsimTNCode: [5,1,3]"

julia> nkd(qs_code)
(5, 1, 3)

julia> validate(qs_code)  # no error indicates success
```
"""
struct QecsimTNCode <: StabilizerCode
    _stabilizers::BitMatrix
    _logical_xs::BitMatrix
    _logical_zs::BitMatrix
    _nkd::Tuple{Int,Int,Union{Int,Missing}}
    _label::String
    tn_code::TensorNetworkCode
end
Model.stabilizers(code::QecsimTNCode) = code._stabilizers
Model.logical_xs(code::QecsimTNCode) = code._logical_xs
Model.logical_zs(code::QecsimTNCode) = code._logical_zs
Model.nkd(code::QecsimTNCode) = code._nkd
Model.label(code::QecsimTNCode) = code._label
function QecsimTNCode(code::TensorNetworkCode; distance=missing, label=nothing)
    # derived defaults
    dnkd = (num_qubits(code), length(code.logicals)÷2, distance)
    dlabel = isnothing(label) ? "QecsimTNCode: [$(dnkd[1]),$(dnkd[2]),$(dnkd[3])]" : label
    return QecsimTNCode(
        _tnpauli_to_bsf(code.stabilizers),
        _tnpauli_to_bsf(code.logicals[1:2:end]),
        _tnpauli_to_bsf(code.logicals[2:2:end]),
        dnkd,
        dlabel,
        code
    )
end

"""
    _tnpauli_to_bsf(pauli::AbstractVector{Int}) -> BitVector
    _tnpauli_to_bsf(paulis::AbstractVector{<:AbstractVector{Int}}) -> BitMatrix

Return the binary symplectic representation of Pauli given in integer vector representation.
"""
function _tnpauli_to_bsf(p::AbstractVector{Int})
    return vcat((x-> x ∈ (1, 2)).(p), (x-> x ∈ (2, 3)).(p))
end
function _tnpauli_to_bsf(ps::AbstractVector{<:AbstractVector{Int}})
    return transpose(reduce(hcat, _tnpauli_to_bsf.(ps)))
end
"""
    _bsf_to_tnpauli(b::AbstractVector{Bool}) -> Vector{Int}
    _bsf_to_tnpauli(bs::AbstractMatrix{Bool}) -> Vector{Vector{Int}}

Return the integer vector representation of Pauli given in binary symplectic representation.
"""
function _bsf_to_tnpauli(b::AbstractVector{Bool})
    return map(x -> (0,1,3,2)[x+1], b[1:end÷2] + 2b[end÷2+1:end])
end
function _bsf_to_tnpauli(bs::AbstractMatrix{Bool})
    return _bsf_to_tnpauli.(b for b in eachrow(bs))
end

"""
    QecsimTNDecoder <: Decoder

    QecsimTNDecoder(chi::Union{Nothing,Integer}=nothing)

Qesim decoder implementation based on the [TNDecode module](@ref). A null value of `chi`
corresponds to exact contraction, otherwise chi defines the bond dimension used in MPS/MPO
contraction.

An `ArgumentError` is thrown if chi is not null or positive. The Qecsim `label` method
returns a label constructed from `chi`. The Qecsim `decode` method expects the keyword
argument `p` for error probability, which defaults to `0.1`, and returns a `DecodeResult`
object with `custom_values` containing the success probability. The success probability is
defined as the ratio of the probability of the mostly-likely logical coset to the sum of the
probabilities of all logical cosets, for example:

    DecodeResult.custom_values = [[0.92]] # e.g. where 0.92 is the success probability

The success probability is added to `custom_values` as a vector so that as Qecsim
`App.qec_run` aggregates success probabilities over many runs it extends the vector rather
than simply summing over values.

!!! note
    Currently the [TNDecode module](@ref) only supports codes that can be laid out on a
    square lattice.

# Examples
```jldoctest
julia> using TensorNetworkCodes.QecsimAdaptors

julia> using Qecsim, Qecsim.GenericModels

julia> using Random:MersenneTwister  # use RNG for reproducible example

julia> code = QecsimTNCode(rotated_surface_code(3); distance=3);

julia> error_model = DepolarizingErrorModel();

julia> decoder = QecsimTNDecoder(4);

julia> label(decoder)
"QecsimTNDecoder (chi=4)"

julia> result = qec_run_once(code, error_model, decoder, 0.1, MersenneTwister(11))
RunResult{Vector{Vector{Float64}}}(true, 1, Bool[0, 0], [[0.9249813981321254]])

julia> success_flag = result.success
true

julia> success_probability = result.custom_values[1][1]
0.9249813981321254
```
"""
struct QecsimTNDecoder <: Decoder
    chi::Union{Nothing,Int}
    function QecsimTNDecoder(chi=nothing)
        if !(isnothing(chi) || chi > 0)
            throw(ArgumentError("QecsimTNDecoder supports `nothing` or positive `chi`."))
        end
        return new(chi)
    end
end
function Model.label(decoder::QecsimTNDecoder)
    chi = decoder.chi
    return isnothing(chi) ? "QecsimTNDecoder" : "QecsimTNDecoder (chi=$(chi))"
end
function Model.decode(
    decoder::QecsimTNDecoder, code::QecsimTNCode, syndrome::AbstractVector{Bool};
    p=0.1, kwargs...
)
    chi = decoder.chi
    tn_syndrome::Vector{Int} = syndrome  # promote to Vector{Int}
    tn_contract_fn = isnothing(chi) ? basic_contract() : mps_contract(chi)
    tn_recovery, predicted_success_rate = tn_decode(
        code.tn_code, tn_syndrome, p;
        contract_fn=tn_contract_fn
    )
    recovery = _tnpauli_to_bsf(tn_recovery)
    # add predicted_success_rate as vector so all values retained (not summed) over runs
    return DecodeResult(recovery=recovery, custom_values=[[predicted_success_rate]])
end

end

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


struct QecsimTNDecoder <: Decoder end

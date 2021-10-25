"""
    find_distance_logicals(code::Quantum_code; max_distance=5) -> Int, Vector{Vector{Int}}

Return the distance of the code and all minimum-weight logical operators.

This method works by brute force. It searches for operators of increasing weight so it works
well for low-distance codes but will be slow for high-distance codes. If during the search
`max_distance` is exceeded then an `ErrorException` is thrown.

# Examples
```jldoctest
julia> d, ls = find_distance_logicals(five_qubit_code());

julia> d, length(ls), ls[1]  # distance, number and example of minimum-weight logicals
(3, 30, [1, 2, 1, 0, 0])
```
"""
function find_distance_logicals(code::QuantumCode; max_distance=5)
    n = num_qubits(code)
    r = length(code.stabilizers)

    lowest_weight_logicals = Vector{Int}[]

    if n == 0  # no qubits so distance is 0
        return 0, lowest_weight_logicals
    end
    if n == r  # no logicals, so distance is n (or ∞?)
        return n, lowest_weight_logicals
    end

    # find distance and lowest_weight_logicals
    distance = n - 1
    locations_iterator = combinations(1:n)
    for locations in locations_iterator

        L = length(locations)

        if L > distance && length(lowest_weight_logicals) > 0  # we have finished
            return distance, lowest_weight_logicals
        end

        if L > max_distance  # we have exceeded max_distance, so throw exception
            error("distance > max_distance = $max_distance; stopping search")
        end

        for l in 1:4^L - 1
            paulis = digits!(zeros(Int, L), l, base=4) # a nice iterator would be better

            operator = zeros(Int, n)
            for index in 1:L
                operator[locations[index]] = paulis[index]
            end

            if all(==(0), pauli_commutation.(Ref(operator), code.stabilizers))
                if (operator in code.logicals
                        || any(==(1), pauli_commutation.(Ref(operator), code.logicals)))
                    distance = L  # = pauli_weight(operator)
                    push!(lowest_weight_logicals, operator)
                end
            end
        end
    end
end

"""
    _find_product_indices(
        operators::AbstractVector{<:AbstractVector{Int}},
        target_operator::AbstractVector{Int}
    ) -> Vector{Int}

Return the indices of the smallest number of operators that as a product give the target
operator. If no such product exists an `ErrorException` is thrown.
"""
function _find_product_indices(operators, target_operator)
    identity = zeros(Int, length(target_operator))
    if target_operator == identity
        return Int[]
    end
    for i in combinations(collect(enumerate(operators)))
        ids, ops = collect(zip(i...))
        op = reduce((a, b)->pauli_product.(a, b), ops; init=identity)
        if op == target_operator
            return collect(ids)
        end
    end
    error("no product of operators yields operator")
end

"""
    find_pure_error(code::QuantumCode, syndrome::AbstractVector{Int}) -> AbstractVector{Int}

Return a pure error which yields the given syndrome with the given code.

The pure error is formed from a product of `code.pure_errors` and is not unique nor
necessarily the lowest-weight error corresponding to the syndrome.

See also [`find_syndrome`](@ref).

# Examples
```jldoctest
julia> code = five_qubit_code();

julia> syndrome = find_syndrome(code, [0, 1, 3, 0, 1])  # error = IXZIX
4-element Vector{Int64}:
 1
 0
 0
 1

julia> pure_error = find_pure_error(code, syndrome)
5-element Vector{Int64}:
 1
 1
 0
 0
 0

julia> find_syndrome(code, pure_error) == syndrome
true
```
"""
function find_pure_error(code::QuantumCode, syndrome::AbstractVector{Int})
    (length(code.pure_errors) == length(syndrome)) || error("invalid syndrome length")
    return pauli_product_pow(code.pure_errors, syndrome)
end

"""
    _find_pure_errors_disordered(stabilizers::AbstractVector{<:AbstractVector{Int}})
        -> Vector{Vector{Int}}

Return pure errors corresponding to a list of stabilizers. The returned list of pure errors
can be used (combined) to find a list of pure errors such that each pure error anticommutes
with precisely one stabilizer; this can be achieved by passing the output of this method to
[`_fix_pure_errors!`](@ref).

See also [`_fix_pure_errors!`](@ref), [`find_pure_errors`](@ref).
"""
function _find_pure_errors_disordered(stabilizers::AbstractVector{<:AbstractVector{Int}})
    stabilizers = copy(stabilizers)  # avoid modifying incoming stabilizers
    pure_errors = Vector{Vector{Int}}()
    r = length(stabilizers)
    if r == 0  # empty stabilizers, so return empty pure_errors
        return pure_errors
    end
    n = length(stabilizers[1])
    remaining = collect(1:r)

    for qubit in 1:r  # TODO: should be 1:n ?
        paulis = [1,2,3]
        indices = Int[]
        for α in remaining
            if stabilizers[α][qubit] in paulis
                setdiff!(paulis, stabilizers[α][qubit])
                push!(indices, α)
            end
            if length(paulis) == 1
                new_pure_error1 = zeros(Int64, n)
                new_pure_error2 = zeros(Int64, n)
                new_pure_error1[qubit] = stabilizers[indices[end]][qubit]
                new_pure_error2[qubit] = stabilizers[indices[end - 1]][qubit]
                push!(pure_errors, new_pure_error1)
                push!(pure_errors, new_pure_error2)
                break
            end
            if α == remaining[end] && length(paulis) == 2
                new_pure_error1 = zeros(Int64, n)
                new_pure_error1[qubit] = paulis[1]
                push!(pure_errors, new_pure_error1)
            end
        end

        remaining = setdiff(remaining, indices)
        for α in remaining
            if stabilizers[α][qubit] == 0
                continue
            end
            new_operator = pauli_product.(stabilizers[α], stabilizers[indices[1]])
            if new_operator[qubit] == 0
                stabilizers[α] = deepcopy(new_operator)
                continue
            end
            if length(paulis) == 2
                continue
            end
            new_operator = pauli_product.(stabilizers[α], stabilizers[indices[2]])
            if new_operator[qubit] == 0
                stabilizers[α] = deepcopy(new_operator)
                continue
            end
            new_operator = pauli_product.(stabilizers[α], stabilizers[indices[1]])
            new_operator = pauli_product.(new_operator, stabilizers[indices[2]])
            if new_operator[qubit] == 0
                stabilizers[α] = deepcopy(new_operator)
                continue
            end
        end
    end

    return pure_errors
end

"""
    _fix_pure_errors!(
        pure_errors::AbstractVector{<:AbstractVector{Int}},
        stabilizers::AbstractVector{<:AbstractVector{Int}}
    )

Reorder and take products of the pure errors so that each pure error anticommutes with
precisely one stabilizer and the order of pure errors respects that of the stabilizers. An
`ErrorException` is thrown if the function cannot succeed; for example if the stabilizers
are not linearly independent.

See also [`_find_pure_errors_disordered`](@ref), [`find_pure_errors`](@ref).
"""
function _fix_pure_errors!(
    pure_errors::AbstractVector{<:AbstractVector{Int}},
    stabilizers::AbstractVector{<:AbstractVector{Int}}
)
    r = length(pure_errors)

    success = 0
    for α in 1:r
        for β in α:r
            if pauli_commutation(pure_errors[β], stabilizers[α]) != 0
                pure_errors[α], pure_errors[β] = pure_errors[β], pure_errors[α]
                success += 1
                break
            end
        end
        for β in 1:r
            if α != β && pauli_commutation(pure_errors[β], stabilizers[α]) != 0
                pure_errors[β] = pauli_product.(pure_errors[α], pure_errors[β])
            end
        end
    end

    # Check if enough pure_errors were found
    if success != r
        error("pure errors cannot be fixed to correspond to stabilizers")
    end
end

"""
    find_pure_errors(stabilizers::AbstractVector{<:AbstractVector{Int}})
        -> Vector{Vector{Int}}

Return pure errors corresponding to a list of stabilizers, such that each pure error
anticommutes with precisely one stabilizer and the order of pure errors respects that of the
stabilizers.

This function is efficient but does not give lowest weight pure errors (you cannot have both
of these properties). An `ErrorException` is thrown if the function cannot succeed; for
example if the stabilizers are not linearly independent.

# Examples
```jldoctest
julia> stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]];

julia> pure_errors = find_pure_errors(stabilizers)
4-element Vector{Vector{Int64}}:
 [0, 1, 0, 0, 0]
 [1, 3, 0, 0, 0]
 [3, 1, 0, 0, 0]
 [1, 0, 0, 0, 0]

julia> [pauli_commutation(s, p) for s in stabilizers, p in pure_errors]  # commutations
4×4 Matrix{Int64}:
 1  0  0  0
 0  1  0  0
 0  0  1  0
 0  0  0  1
```
"""
function find_pure_errors(stabilizers::Array{Array{Int64,1},1})
    pure_errors = _find_pure_errors_disordered(stabilizers)
    _fix_pure_errors!(pure_errors, stabilizers)
    return pure_errors
end

"""
    find_syndrome(code::QuantumCode, error::AbstractVector{Int}) -> AbstractVector{Int}

Return the syndrome yielded by the given error with the given code.

The syndrome is a list of 1 and 0 of the same length as `code.stabilizers`, where 1
indicates the error anticommutes with the corresponding stabilizer.

See also [`find_pure_error`](@ref).

# Examples
```jldoctest
julia> syndrome = find_syndrome(five_qubit_code(), [0, 1, 3, 0, 1])  # error = IXZIX
4-element Vector{Int64}:
 1
 0
 0
 1
```
"""
function find_syndrome(code::QuantumCode, error_operator::AbstractVector{Int})
    num_qubits(code) == length(error_operator) || error("invalid error operator length")
    return pauli_commutation.(code.stabilizers, Ref(error_operator))
end

"""
    gauge_code(code::SimpleCode, logical_qubit::Int, logical_pauli::Int) -> SimpleCode

Given a simple code with ``k`` logicals on ``n`` physical qubits, return a new simple code
with ``k - 1`` logicals on ``n`` physical qubits by adding a logical operator as a
stabilizer, where `logical_qubit` indexes which logical qubit is gauged and `logical_pauli`
indicates which logical Pauli is added to the stabilizers.

A `ErrorException` is thrown if `logical_qubit` indexes a non-existant logical qubit, or if
`logical_pauli` is not in `1:3` (logical identity does not fix a gauge).

# Examples
```jldoctest
julia> code = five_qubit_code();

julia> code.stabilizers
4-element Vector{Vector{Int64}}:
 [1, 3, 3, 1, 0]
 [0, 1, 3, 3, 1]
 [1, 0, 1, 3, 3]
 [3, 1, 0, 1, 3]

julia> code.logicals
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1]
 [3, 3, 3, 3, 3]

julia> new_code = gauge_code(code, 1, 3);  # gauge logical qubit 1 using logical Z

julia> new_code.stabilizers
5-element Vector{Vector{Int64}}:
 [1, 3, 3, 1, 0]
 [0, 1, 3, 3, 1]
 [1, 0, 1, 3, 3]
 [3, 1, 0, 1, 3]
 [3, 3, 3, 3, 3]

julia> new_code.logicals
Vector{Int64}[]
```
"""
function gauge_code(code::SimpleCode, logical_qubit::Int, logical_pauli::Int)
    k = Int(length(code.logicals) / 2)
    # preconditions
    (logical_qubit in 1:k) || error("logical qubit index out of bounds!")
    (logical_pauli in 0:3) || error("unknown logical pauli operator!")
    (logical_pauli != 0) || error("logical identity doesn't fix a gauge!")
    # new code fields
    name = "$logical_qubit/$(pauli_rep_change(logical_pauli)) gauged $(code.name)"
    output_stabilizers = deepcopy(code.stabilizers)
    output_logicals = deepcopy(code.logicals)
    output_pure_errors = deepcopy(code.pure_errors)
    # remove gauged logicals
    gauged_logicals = [
        popat!(output_logicals, 2*logical_qubit-1),
        popat!(output_logicals, 2*logical_qubit-1),
    ]
    # new stabilizer from gauged logicals
    logical_powers = [[0, 0], [1, 0], [1, 1], [0, 1]][logical_pauli + 1]
    new_stabilizer = pauli_product_pow(gauged_logicals, logical_powers)
    push!(output_stabilizers, new_stabilizer)
    # new pure error from gauged logicals
    logical_powers = [[0, 0], [1, 1], [0, 1], [1, 1]][logical_pauli + 1]
    new_pure_error = pauli_product_pow(gauged_logicals, logical_powers)
    push!(output_pure_errors, new_pure_error)
    _fix_pure_errors!(output_pure_errors, output_stabilizers)

    return SimpleCode(name, output_stabilizers, output_logicals, output_pure_errors)
end

"""
    num_qubits(code::QuantumCode) -> Int

Return the number of physical qubits of the code.
"""
function num_qubits(code::QuantumCode)
    return length(code.stabilizers) == 0 ? 0 : length(code.stabilizers[1])
end

"""
    permute_code(code::SimpleCode, permutation) -> SimpleCode

Return a new simple code with the physical qubits permuted relative to the given code,
according to the permutation.

The `permutation` is expected in the format used for `Base.permute!` and it is applied to
each stabilizer, logical and pure error of the code. No checking is done to verify that
`permuation` is a valid.

# Examples
```jldoctest
julia> code = five_qubit_code();

julia> code.name
"Five qubit code"

julia> code.stabilizers
4-element Vector{Vector{Int64}}:
 [1, 3, 3, 1, 0]
 [0, 1, 3, 3, 1]
 [1, 0, 1, 3, 3]
 [3, 1, 0, 1, 3]

julia> new_code = permute_code(code, [2, 1, 3, 4, 5]);

julia> new_code.name
"Five qubit code [2, 1, 3, 4, 5]"

julia> new_code.stabilizers
4-element Vector{Vector{Int64}}:
 [3, 1, 3, 1, 0]
 [1, 0, 3, 3, 1]
 [0, 1, 1, 3, 3]
 [1, 3, 0, 1, 3]
```
"""
function permute_code(code::SimpleCode, permutation)
    output_code = SimpleCode("$(code.name) $(permutation)",
        deepcopy(code.stabilizers),
        deepcopy(code.logicals),
        deepcopy(code.pure_errors)
    )
    permute!.(output_code.stabilizers, Ref(permutation))
    permute!.(output_code.logicals, Ref(permutation))
    permute!.(output_code.pure_errors, Ref(permutation))
    return output_code
end

"""
    purify_code(code::SimpleCode) -> SimpleCode

Given a simple code with ``k`` logicals on ``n`` physical qubits, return a new simple code
with ``0`` logicals on ``n + k`` physical qubits.

# Examples
```jldoctest
julia> code = purify_code(five_qubit_code());

julia> num_qubits(code), length(code.logicals)
(6, 0)
```
"""
function purify_code(code::SimpleCode)
    g = code.stabilizers
    l = code.logicals
    K = length(l)
    k = Int(K / 2)
    n = num_qubits(code)

    output_stabilizers = vcat.(Ref(zeros(Int, k)), g)
    new_stabilizers = vcat.(Ref(zeros(Int, k)), l)
    output_pure_errors = vcat.(Ref(zeros(Int, k)), code.pure_errors)
    new_pure_errors = [zeros(Int, n + k) for _ in 1:K]

    for α in 1:length(new_stabilizers)
        if isodd(α)
            β = Int((α + 1) ÷ 2)
            new_stabilizers[α][β] = 1
            new_pure_errors[α][β] = 3
        elseif iseven(α)
            β = Int(α ÷ 2)
            new_stabilizers[α][β] = 3
            new_pure_errors[α][β] = 1
        end
    end

    name = "Purified $(code.name)"
    output_stabilizers = vcat(output_stabilizers, new_stabilizers)
    output_logicals = Vector{Int}[]
    output_pure_errors = vcat(output_pure_errors, new_pure_errors)
    _fix_pure_errors!(output_pure_errors, output_stabilizers)

    return SimpleCode(name, output_stabilizers, output_logicals, output_pure_errors)
end

"""
    verify_code(code::QuantumCode; log_warn=true) -> Bool

Return true if the code satisfied the properties of a valid code, or false otherwise. if the
code is not valid and `log_warn` is true then a warning is logged with the specific reason.

The following checks are performed:
* Number of stabilizers, pure errors and logicals are consistent.
* Stabilizers are independent and mutually commute.
* Pure errors anticommute with corresponding stabilizers and commute with other stabilizers.
* Logicals commute with stabilizers.

The following checks are not yet performed:
* Logical commutation relations.
"""
function verify_code(code::QuantumCode; log_warn=true)
    n = num_qubits(code)
    r = length(code.stabilizers)
    p = length(code.pure_errors)
    l = length(code.logicals)

    # Check we have the right number of operators.
    if p != r
        log_warn && @warn "number of stabilizers and pure errors don't match!"
        return false
    end
    if n != r + l / 2
        log_warn && @warn "number of stabilizers and logicals don't add up!"
        return false
    end
    # Check stabilizers are independent.
    if !pauli_are_independent(code.stabilizers)
        log_warn && @warn "stabilizers aren't independent!"
        return false
    end
    # Check stabilizers mutually commute.
    if !pauli_are_commuting(code.stabilizers)
        log_warn && @warn "stabilizers don't commute!"
        return false
    end
    # Check pure errors anticommute with corresponding stabilizers.
    if !all(==(1), pauli_commutation.(code.stabilizers, code.pure_errors))
        log_warn && @warn "pure errors don't anticommute with corresponding stabilizers!"
        return false
    end
    # Check pure error do not anticommute with other stabilizers.
    for α in 1:r, β in 1:r
        if α != β && pauli_commutation(code.stabilizers[α], code.pure_errors[β]) == 1
            log_warn && @warn "pure errors anticommute with the wrong stabilizers!"
            return false
        end
    end
    # Check logicals commute with stabilizers.
    for logical in code.logicals
        if !all(==(0), pauli_commutation.(code.stabilizers, Ref(logical)))
            log_warn && @warn "logicals don't commute with stabilizers!"
            return false
        end
    end
    # TODO: check logical commutation relations.

    return true
end

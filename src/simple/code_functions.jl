"""
    distance_logicals(code::Quantum_code; max_distance=5) -> Int, Vector{Vector{Int}}

Return the distance of the code and all minimum-weight logical operators.

This method works by brute force. It searches for operators of increasing weight so it works
well for low-distance codes but will be slow for high-distance codes. If during the search
`max_distance` is exceeded then an `ErrorException` is thrown.

# Examples
```jldoctest
julia> d, ls = distance_logicals(five_qubit_code());

julia> d, length(ls), ls[1]  # distance, number and example of minimum-weight logicals
(3, 30, [1, 2, 1, 0, 0])
```
"""
function distance_logicals(code::QuantumCode; max_distance=5)
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
    num_qubits(code::QuantumCode) -> Int

Return the number of physical qubits of the code.
"""
function num_qubits(code::QuantumCode)
    return length(code.stabilizers) == 0 ? 0 : length(code.stabilizers[1])
end

"""
    verify_code(code::QuantumCode, log_warn=true) -> Bool

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
function verify_code(code::QuantumCode, log_warn=true)
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

"""
    permute(code,permutation)

Permutes physical qubits of a `SimpleCode` returning a new `SimpleCode`.
"""
function permute(
        code::SimpleCode,
        permutation::Array{Int64})

    output_code =
    SimpleCode(code.name * " " * string(permutation),
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
    fix_pure_errors!(pure_errors,stabilizers)

After modifying stabilizers (e.g, taking products to get a new generating
set), pure_errors may not satisfy the correct anticommutation relations.
This fixes this, returning true if it succeeded and false if it was not
possible.
"""
function fix_pure_errors!(
        pure_errors::Array{Array{Int64,1},1},
        stabilizers::Array{Array{Int64,1},1})

    success = 0
    r = length(pure_errors)


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
        return false
    else
        return true
    end
end





"""
    generate_disordered_pure_errors(stabilizers)

This finds pure_errors corresponding to a set of stabilizers.  The
pure errors don't necessarily have the same order as the stabilizers.
"""
function generate_disordered_pure_errors(stabilizers::Array{Array{Int64,1},1})

    r = length(stabilizers)
    n = length(stabilizers[1])
    remaining = collect(1:r)
    pure_errors = Array{Array{Int64,1},1}()

    for qubit in 1:r
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
    generate_pure_errors(stabilizers)

This finds pure_errors corresponding to a set of stabilizers.  It is
efficient but does not give lowest weight pure errors (you actually
can't have both of these properties).
"""
function generate_pure_errors(stabilizers::Array{Array{Int64,1},1})

    output = generate_disordered_pure_errors(stabilizers)
    fix_pure_errors!(output, stabilizers)

    return output
end





















"""
    low_weight_stabilizers(stabilizers)

Finds a new set of stabilizers that have low weight by first
finding the lowest weight element of the group, then the second
lowest independent element and so on.  This works by brute force.

Warning: this doesn't update/fix `pure_errors`!
"""
function low_weight_stabilizers(stabilizers::Array{Array{Int64,1}})

    g = stabilizers
    n = length(g[1])
    r = length(g)

    new_stabilizers = Array{Int64,1}[]

    for α in 1:r
        op_weight = n + 1
        lowest_weight_stabilizer = zeros(Int, n)
        for m in 1:2^r - 1
            powers = digits!(zeros(Int8, r), m, base=2)

            operator = pauli_product_pow(stabilizers, powers)
            w = pauli_weight(operator)
            temp = deepcopy(new_stabilizers)

            if (w < op_weight && pauli_are_independent(push!(temp, operator))
                    && operator != zeros(Int, n))
                op_weight = w
                lowest_weight_stabilizer = operator
            end

            if op_weight <= 4
                @goto find_next
            end
        end

        @label find_next
        push!(new_stabilizers, lowest_weight_stabilizer)
    end

    return new_stabilizers
end



low_weight_stabilizers(code::QuantumCode) = low_weight_stabilizers(code.stabilizers)





"""
    are_physically_equivalent(x::QuantumCode,y::QuantumCode)

Checks if two `QuantumCodes` are equal in a physical sense, i.e., do all
stabilizers commute, do logicals on `x` commute with stabilizers of `y`, etc.
"""
function are_physically_equivalent(x::QuantumCode, y::QuantumCode)
    if num_qubits(x) != num_qubits(y)
        return false
    end
    if length(x.stabilizers) != length(y.stabilizers)
        return false
    end
    if length(x.logicals) != length(y.logicals)
        return false
    end


    if !pauli_are_commuting(vcat(x.stabilizers, y.stabilizers))
        return false
    end
    for logical in x.logicals
        if !pauli_are_commuting(vcat(y.stabilizers, [logical]))
            return false
        end
    end
    for logical in y.logicals
        if !pauli_are_commuting(vcat(x.stabilizers, [logical]))
            return false
        end
    end

    # Check that the stabilizer groups are equal
    for stabilizer in y.stabilizers
        if pauli_are_independent(vcat(x.stabilizers, [stabilizer]))
            return false
        end
    end


    # Check if each logical of x anticommutes with only one logical
    # of y, and no two logicals of x anticommute with the same logical of
    # x.
    new_logicals = deepcopy(x.logicals)
    if !fix_pure_errors!(new_logicals, y.logicals)
        return false
    end

    return true
end





"""
    purify_code(code)

Given a `SimpleCode` with `k` logicals on `n` physical qubits, returns
a `SimpleCode` with 0 logicals on `n+k` physical qubits.
"""
function purify_code(code::SimpleCode)
    g = code.stabilizers
    l = code.logicals
    K = length(l)
    k = Int(K / 2)
    n = num_qubits(code)


    output_stabilizers = vcat.(Ref(zeros(Int64, k)), g)
    new_stabilizers = vcat.(Ref(zeros(Int64, k)), l)
    output_pure_errors = vcat.(Ref(zeros(Int64, k)), code.pure_errors)
    new_pure_errors = [zeros(Int64, n + k) for _ in 1:K]


    for α in 1:length(new_stabilizers)
        if isodd(α)
            β = Int((α + 1) / 2)
            new_stabilizers[α][β] = 1
            new_pure_errors[α][β] = 3
        elseif iseven(α)
            β = Int(α / 2)
            new_stabilizers[α][β] = 3
            new_pure_errors[α][β] = 1
        end
    end


    output_stabilizers = vcat(output_stabilizers, new_stabilizers)
    output_logicals = []
    output_pure_errors = vcat(output_pure_errors, new_pure_errors)
    fix_pure_errors!(output_pure_errors, output_stabilizers)

    name = "Pureified " * code.name

    return SimpleCode(name, output_stabilizers, output_logicals, output_pure_errors)
end





"""
    gauge_code(code,logical_power_list,which_logicals)

Given a `SimpleCode` with `k` logicals on `n` physical qubits, returns
a `SimpleCode` with `k - length(logical_power_list)/2` logicals on `n`
physical qubits by adding logical operators as stabilizer.
"""
function gauge_code(
        code::SimpleCode,
        logical_power_list::Array{Array{Int64,1},1},
        which_logicals::Array{Int64,1})

    l = code.logicals
    K = length(l)

    if length(which_logicals) != length(logical_power_list)
        error("incorrect number of powers!")
    end
    if [0,0] ∈ logical_power_list
        error("this doesn't fix a gauge!")
    end


    output_stabilizers = deepcopy(code.stabilizers)
    output_pure_errors = deepcopy(code.pure_errors)

    for β in 1:length(which_logicals)
        α = which_logicals[β]
        logical_powers = copy(logical_power_list[β])
        op_range = 2 * (α - 1) + 1:2 * (α)

        new_stabilizer = pauli_product_pow(l[op_range], logical_powers)
        if logical_powers[2] == 1
            logical_powers[1] = mod(logical_powers[1] + 1, 2)
        elseif logical_powers[1] == 1
            logical_powers[2] = mod(logical_powers[2] + 1, 2)
        end
        new_pure_error = pauli_product_pow(l[op_range], logical_powers)

        push!(output_stabilizers, new_stabilizer)
        push!(output_pure_errors, new_pure_error)
    end


    output_logicals = Array{Int64,1}[]
    for α in 1:Int(K / 2)
        if !(α ∈ which_logicals)
            push!(output_logicals, code.logicals[2 * (α - 1) + 1])
            push!(output_logicals, code.logicals[2 * (α)])
        end
    end


    fix_pure_errors!(output_pure_errors, output_stabilizers)
    name = string(logical_power_list) * " gauged " * code.name

    return SimpleCode(name, output_stabilizers, output_logicals, output_pure_errors)
end





"""
    random_stabilizer_state(n)

Returns a random `SimpleCode` with `0` logicals and `n` physical
qubits.
"""
function random_stabilizer_state(n::Int64)

    stabilizers = Array{Int64,1}[]
    for α in 1:n
        @label try_again
        new_stabilizer = rand(0:3, n)

        if pauli_weight(new_stabilizer) <= 1
            @goto try_again
        elseif !pauli_are_independent(vcat(stabilizers, [new_stabilizer]))
            @goto try_again
        elseif !pauli_are_commuting(vcat(stabilizers, [new_stabilizer]))
            @goto try_again
        end

        push!(stabilizers, new_stabilizer)
    end


    # Next find pure_errors
    pure_errors = generate_pure_errors(stabilizers)
    logicals = []

    return SimpleCode(
        "random stabilizer state",
        stabilizers,
        logicals,
        pure_errors)
end





"""
    random_code(n,k)

Returns a random `SimpleCode` with `k` logicals and `n` physical
qubits.
"""
function random_code(n::Int64, k::Int64)

    if k > n
        error("k > n!")
    end

    r = n - k

    stabilizer_state = random_stabilizer_state(n)
    stabilizers = stabilizer_state.stabilizers[1:r]
    pure_errors = stabilizer_state.pure_errors[1:r]

    logicals = Array{Int64,1}[]
    for α in r + 1:n
        push!(logicals, stabilizer_state.stabilizers[α])
        push!(logicals, stabilizer_state.pure_errors[α])
    end


    return SimpleCode(
        "random code",
        stabilizers,
        logicals,
        pure_errors)
end

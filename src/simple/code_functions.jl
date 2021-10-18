"""
    size(code)

Returns the size (number of physical qubits) of a `Quantum_code`.
"""
function Base.size(code::QuantumCode)
    if length(code.stabilizers) == 0
        return 0
    end

    return length(code.stabilizers[1])
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

    permute!.(output_code.stabilizers,Ref(permutation))
    permute!.(output_code.logicals,Ref(permutation))
    permute!.(output_code.pure_errors,Ref(permutation))

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
            if pauli_commutation(pure_errors[β],stabilizers[α]) != 0
                pure_errors[α],pure_errors[β] = pure_errors[β],pure_errors[α]
                success += 1
                break
            end
        end
        for β in 1:r
            if α != β && pauli_commutation(pure_errors[β],stabilizers[α]) != 0
                pure_errors[β] = pauli_product.(pure_errors[α],pure_errors[β])
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
                setdiff!(paulis,stabilizers[α][qubit])
                push!(indices,α)
            end
            if length(paulis) == 1

                new_pure_error1 = zeros(Int64,n)
                new_pure_error2 = zeros(Int64,n)
                new_pure_error1[qubit] = stabilizers[indices[end]][qubit]
                new_pure_error2[qubit] = stabilizers[indices[end-1]][qubit]
                push!(pure_errors,new_pure_error1)
                push!(pure_errors,new_pure_error2)
                break
            end
            if α == remaining[end] && length(paulis) == 2

                new_pure_error1 = zeros(Int64,n)
                new_pure_error1[qubit] = paulis[1]
                push!(pure_errors,new_pure_error1)
            end
        end


        remaining = setdiff(remaining,indices)
        for α in remaining
            if stabilizers[α][qubit] == 0
                continue
            end
            new_operator = pauli_product.(stabilizers[α],stabilizers[indices[1]])
            if new_operator[qubit] == 0
                stabilizers[α] = deepcopy(new_operator)
                continue
            end
            if length(paulis) == 2
                continue
            end
            new_operator = pauli_product.(stabilizers[α],stabilizers[indices[2]])
            if new_operator[qubit] == 0
                stabilizers[α] = deepcopy(new_operator)
                continue
            end
            new_operator = pauli_product.(stabilizers[α],stabilizers[indices[1]])
            new_operator = pauli_product.(new_operator,stabilizers[indices[2]])
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
    fix_pure_errors!(output,stabilizers)

    return output
end





"""
    verify_code(code::QuantumCode)

Checks most properties of a quantum error correcting code to make sure it's
sensible, e.g., do the stabilizers commute.
"""
function verify_code(code::QuantumCode)

    n = size(code)
    r = length(code.stabilizers)
    p = length(code.pure_errors)
    l = length(code.logicals)


    # Do we have the right number of operators?
    if p != r
        println("number of stabilizers and pure errors do not match!")
        return false
    end
    if n != r + l/2
        println("number of operators do not add up!")
        return false
    end


    # Do pure errors fulfil their role?
    if pauli_commutation.(code.stabilizers,code.pure_errors) != ones(Int64,r)
        println("pure errors don't anticommute with their corresponding
            stabilizers!")
        return false
    end

    for α in 1:r, β in 1:r
        if α != β && pauli_commutation(code.stabilizers[α],code.pure_errors[β]) == 1
            println("pure errors anticommute with the wrong stabilizers!")
            return false
        end
    end


    # Do stabilizers commute?
    if do_they_commute(code.stabilizers) != 0
        println("stabilizers don't commute!")
        return false
    end


    # Are stabilizers independent?  (This may be redundant)
    if !are_they_independent(code.stabilizers)
        println("stabilizers aren't independent!")
        return false
    end


    # Do logicals commute with stabilizers?
    for logical in code.logicals
        operators = deepcopy(code.stabilizers)
        push!(operators,logical)
        if do_they_commute(operators) != 0
            println("logicals don't commute with stabilizers!")
            return false
        end
    end

    # Should check logical commutation relations also.


    return true
end





"""
    distance(code::Quantum_code) -> Int64

Find the distance of a code, i.e., the lowest weight of a nontrivial logical operator.
This works by brute force, but it starts with low-weight operators, so it's good for
low-distance codes.
"""
function distance(
        logicals::Array{Array{Int64,1}},
        stabilizers::Array{Array{Int64,1}};
        max_distance = 5)

    kk = length(logicals)
    n = length(stabilizers[1])
    r = length(stabilizers)


    if r == n
        return n, [Int[]]  # no logicals, so distance is n (or ∞?)
    end


    distance = n - 1
    lowest_weight_logicals = Array{Int64,1}[]
    locations_iterator = combinations(1:n)

    for locations in locations_iterator
        L = length(locations)


        if L > distance && length(lowest_weight_logicals) > 0
            return distance, lowest_weight_logicals
        end
        if L > max_distance
            println("Distance > " * string(max_distance)*".  Stopping now!")
            return nothing
        end


        for l in 1:4^L-1
            paulis = digits!(zeros(Int64,L),l,base = 4) # a nice iterator would be better

            operator = zeros(Int64,n)
            for index in 1:L
                operator[locations[index]] = paulis[index]
            end


            if pauli_commutation.(Ref(operator),stabilizers) == zeros(Int64,r)
                if (pauli_commutation.(Ref(operator),logicals) != zeros(Int64,kk) ||
                    operator ∈ logicals )
                distance = L  # = weight(operator)
                push!(lowest_weight_logicals,operator)
                end
            end
        end
    end
end



function distance(code::QuantumCode;max_distance =5)
    if size(code) == 0
        return 0
    end

    return distance(code.logicals,code.stabilizers;max_distance)
end





# Second method when you don't know logicals
function distance(
        stabilizers::Array{Array{Int64,1}};
        max_distance = 5)

    n = length(stabilizers[1])
    r = length(stabilizers)


    if r == n
        return n, [Int[]]  # no logicals, so distance is n (or ∞?)
    end


    distance = n
    lowest_weight_logicals = Array{Int64,1}[]
    locations_iterator = combinations(1:n)

    for locations in locations_iterator
        L = length(locations)

        if L > max_distance
            return "distance > " * string(max_distance)
        end

        if L > distance && length(lowest_weight_logicals) > 0
            return distance, lowest_weight_logicals
        end

        for l in 1:4^L-1
            paulis = digits!(zeros(Int64,L),l,base = 4) # a nice iterator would be better

            operator = zeros(Int64,n)
            for index in 1:L
                operator[locations[index]] = paulis[index]
            end

            if (pauli_commutation.(Ref(operator),stabilizers) == zeros(Int64,r) &&
                    are_they_independent(vcat(stabilizers,[operator])))
                distance = L
                push!(lowest_weight_logicals,operator)
            end
        end
    end
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
        lowest_weight_stabilizer = zeros(Int,n)
        for m in 1:2^r-1
            powers = digits!(zeros(Int8,r),m,base = 2)

            operator = product_with_powers(stabilizers,powers)
            w = weight(operator)
            temp = deepcopy(new_stabilizers)

            if (w < op_weight && are_they_independent(push!(temp,operator))
                    && operator != zeros(Int,n))
                op_weight = w
                lowest_weight_stabilizer = operator
            end

            if op_weight <= 4
                @goto find_next
            end
        end

        @label find_next
        push!(new_stabilizers,lowest_weight_stabilizer)
    end

    return new_stabilizers
end



low_weight_stabilizers(code::QuantumCode) = low_weight_stabilizers(code.stabilizers)





"""
    are_physically_equivalent(x::QuantumCode,y::QuantumCode)

Checks if two `QuantumCodes` are equal in a physical sense, i.e., do all
stabilizers commute, do logicals on `x` commute with stabilizers of `y`, etc.
"""
function are_physically_equivalent(x::QuantumCode,y::QuantumCode)
    if size(x) != size(y)
        return false
    end
    if length(x.stabilizers) != length(y.stabilizers)
        return false
    end
    if length(x.logicals) != length(y.logicals)
        return false
    end


    if do_they_commute(vcat(x.stabilizers,y.stabilizers)) != 0
        return false
    end
    for logical in x.logicals
        if do_they_commute(vcat(y.stabilizers,[logical])) != 0
            return false
        end
    end
    for logical in y.logicals
        if do_they_commute(vcat(x.stabilizers,[logical])) != 0
            return false
        end
    end

    # Check that the stabilizer groups are equal
    for stabilizer in y.stabilizers
        if are_they_independent(vcat(x.stabilizers,[stabilizer]))
            return false
        end
    end


    # Check if each logical of x anticommutes with only one logical
    # of y, and no two logicals of x anticommute with the same logical of
    # x.
    new_logicals = deepcopy(x.logicals)
    if !fix_pure_errors!(new_logicals,y.logicals)
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
    k = Int(K/2)
    n = size(code)


    output_stabilizers = vcat.(Ref(zeros(Int64,k)),g)
    new_stabilizers = vcat.(Ref(zeros(Int64,k)),l)
    output_pure_errors = vcat.(Ref(zeros(Int64,k)),code.pure_errors)
    new_pure_errors = [zeros(Int64,n+k) for _ in 1:K]


    for α in 1:length(new_stabilizers)
        if isodd(α)
            β = Int((α+1)/2)
            new_stabilizers[α][β] = 1
            new_pure_errors[α][β] = 3
        elseif iseven(α)
            β = Int(α/2)
            new_stabilizers[α][β] = 3
            new_pure_errors[α][β] = 1
        end
    end


    output_stabilizers = vcat(output_stabilizers,new_stabilizers)
    output_logicals = []
    output_pure_errors = vcat(output_pure_errors,new_pure_errors)
    fix_pure_errors!(output_pure_errors,output_stabilizers)

    name = "Pureified " * code.name

    return SimpleCode(name,output_stabilizers,output_logicals,output_pure_errors)
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
        op_range = 2*(α-1)+1:2*(α)

        new_stabilizer = product_with_powers(l[op_range],logical_powers)
        if logical_powers[2] == 1
            logical_powers[1] = mod(logical_powers[1]+1,2)
        elseif logical_powers[1] == 1
            logical_powers[2] = mod(logical_powers[2]+1,2)
        end
        new_pure_error = product_with_powers(l[op_range],logical_powers)

        push!(output_stabilizers,new_stabilizer)
        push!(output_pure_errors,new_pure_error)
    end


    output_logicals = Array{Int64,1}[]
    for α in 1:Int(K/2)
        if !(α ∈ which_logicals)
            push!(output_logicals,code.logicals[2*(α-1)+1])
            push!(output_logicals,code.logicals[2*(α)])
        end
    end


    fix_pure_errors!(output_pure_errors,output_stabilizers)
    name = string(logical_power_list) * " gauged " * code.name

    return SimpleCode(name,output_stabilizers,output_logicals,output_pure_errors)
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
        new_stabilizer = rand(0:3,n)

        if weight(new_stabilizer) <= 1
            @goto try_again
        elseif !are_they_independent(vcat(stabilizers,[new_stabilizer]))
            @goto try_again
        elseif do_they_commute(vcat(stabilizers,[new_stabilizer])) == 1
            @goto try_again
        end

        push!(stabilizers,new_stabilizer)
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
function random_code(n::Int64,k::Int64)

    if k > n
        error("k > n!")
    end

    r = n-k

    stabilizer_state = random_stabilizer_state(n)
    stabilizers = stabilizer_state.stabilizers[1:r]
    pure_errors = stabilizer_state.pure_errors[1:r]

    logicals = Array{Int64,1}[]
    for α in r+1:n
        push!(logicals,stabilizer_state.stabilizers[α])
        push!(logicals,stabilizer_state.pure_errors[α])
    end


    return SimpleCode(
        "random code",
        stabilizers,
        logicals,
        pure_errors)
end

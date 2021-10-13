# import Base: size
# using Combinatorics





"""
    size(code)

Returns the size (number of physical qubits) of a `Quantum_code`.
"""
function size(code::QuantumCode)
    if length(code.stabilizers) == 0
        return 0
    end

    return length(code.stabilizers[1])
end





"""
    pauli_product(a,b) -> Int64

Gives the result of multiplying two Paulis.  Note 0 is the identity (not 4 like
in the holo code notebooks).
"""
function pauli_product(a::Int,b::Int)
    if a == 0
        return b
    elseif b == 0
        return a
    elseif a == b
        return 0
    else
        for n in 1:3
            if n != a && n != b
                return n
            end
        end
    end
end





"""
    do_they_commute(a::Int,b::Int) -> 0 or 1

Checks if two single Paulis commute (0) or don't commute (1).
"""
function do_they_commute(a::Int,b::Int)
    if a == 0
        return 0
    elseif b == 0
        return 0
    elseif a == b
        return 0
    else
        return 1
    end
end



"""
Second method: checks if two Pauli vectors commute (0) or don't commute (1).
"""
function do_they_commute(a::Array{Int64,1},b::Array{Int64,1})
    output = 0
    for n in 1:length(a)
        output += do_they_commute(a[n],b[n])
    end
    return mod(output,2)
end



"""
Third method: checks if a set of Pauli vectors is commuting (0) or not commuting (1).
"""
function do_they_commute(a::Array{Array{Int64,1},1})
    num_operators = length(a)
    for n in 1:num_operators, m in 1:n-1
        if do_they_commute(a[n],a[m]) == 1
            return 1
            break
        end
    end
    return 0
end





"""
    pauli_rep_change(pauli::Char) -> Int

Converts from alphabetical to numerical representation of Paulis.
"""
function pauli_rep_change(pauli::Char)
    if pauli ==  'I'
        return 0
    elseif pauli == 'x'
        return 1
    elseif pauli == 'y'
        return 2
    elseif pauli == 'z'
        return 3
    end
end



"""
    pauli_rep_change(pauli::Int) -> Char

Converts from numerical to alphabetical representation of Paulis.
"""
function pauli_rep_change(pauli::Int)
    if pauli == 0
        return 'I'
    elseif pauli == 1
        return 'x'
    elseif pauli == 2
        return 'y'
    elseif pauli == 3
        return 'z'
    end
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
    product_with_powers(operators,powers) -> Array{Int64,1}

Give a list of operators and binary powers, this returns their
product with the appropriate powers.
"""
function product_with_powers(operators,powers)::Array{Int64,1}
    output = fill(0,length(operators[1]))
    for n in 1:length(operators)
        if powers[n] != 0
            if !(powers[n] in [0,1])
               error("powers of Paulis are not in [0,1]!")
            end
            output = pauli_product.(output,operators[n])
        end
    end
    return output
end

product_with_powers(g,l,powers) = product_with_powers(vcat(g,l),powers)





"""
Find the weight of a vector of Pauli operators, i.e., the number of entries with
non identity elements.
"""
function weight(operator)
    output = 0
    for n in 1:length(operator)
        if operator[n] != 0
            output += 1
        end
    end
    return output
end

"""
    pauli_are_commuting(operators) -> Bool

Return true if the Pauli operators mutually commute, or false otherwise, where `operators`
is an iterable of `AbstractVector{Int}`.

# Examples
```jldoctest
julia> stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]];

julia> pauli_are_commuting(stabilizers)  # stabilizers mutually commute
true

julia> logical_x, logical_z = [1, 1, 1, 1, 1], [3, 3, 3, 3, 3];  # XXXXX, ZZZZZ

julia> pauli_are_commuting([logical_x, logical_z])  # logicals do not commute
false
```
"""
function pauli_are_commuting(operators)
    return all(==(0), pauli_commutation(a, b) for (a, b) in combinations(operators, 2))
end

"""
    pauli_are_independent(operators::AbstractVector{<:AbstractVector{Int}}) -> Bool

Return true if the Pauli operators are linearly independent, or false otherwise.

# Examples
```jldoctest
julia> pauli_are_independent([[3, 3, 0], [0, 3, 3]])  # ZZI, IZZ
true

julia> pauli_are_independent([[3, 3, 0], [0, 3, 3], [3, 0, 3]])  # ZZI, IZZ, ZIZ
false
```
"""
function pauli_are_independent(operators::AbstractVector{<:AbstractVector{Int}})

    num_operators = length(operators)
    num_qubits = length(operators[1])
    remaining = collect(1:num_operators)

    for qubit in 1:num_qubits
        paulis = [1,2,3]
        indices = Int[]
        for α in remaining
            if operators[α][qubit] in paulis
                setdiff!(paulis, operators[α][qubit])
                push!(indices, α)
            end
            if length(paulis) == 1
                break
            end
        end

        remaining = setdiff(remaining, indices)
        for α in remaining
            if operators[α][qubit] == 0
                continue
            end
            new_operator = pauli_product.(operators[α], operators[indices[1]])
            if new_operator[qubit] == 0
                operators[α] = deepcopy(new_operator)
                continue
            end
            if length(paulis) == 2
                continue
            end
            new_operator = pauli_product.(operators[α], operators[indices[2]])
            if new_operator[qubit] == 0
                operators[α] = deepcopy(new_operator)
                continue
            end
            new_operator = pauli_product.(operators[α], operators[indices[1]])
            new_operator = pauli_product.(new_operator, operators[indices[2]])
            if new_operator[qubit] == 0
                operators[α] = deepcopy(new_operator)
                continue
            end
        end

        for α in remaining
            if operators[α] == zeros(Int64, num_qubits)
                return false
            end
        end
    end

    return true

end

"""
    pauli_commutation(a::Int, b::Int) -> 0 or 1
    pauli_commutation(a::AbstractVector{Int}, b::AbstractVector{Int}) -> 0 or 1

Return the commutation relation of two Paulis, i.e. 0 if they commute and 1 if not.

# Examples
```jldoctest
julia> pauli_commutation(1, 1)  # X commutes with itself
0

julia> pauli_commutation(1, 3)  # X does not commute with Z
1
```

```jldoctest
julia> pauli_commutation([1, 3, 3, 1, 0], [1, 1, 1, 1, 1])  # XZZXI commutes with XXXXX
0

julia> stabilizers = [[1, 3, 3, 1, 0], [0, 1, 3, 3, 1], [1, 0, 1, 3, 3], [3, 1, 0, 1, 3]];

julia> error = [1, 1, 0, 0, 0];  # XXIII

julia> syndrome = pauli_commutation.(stabilizers, Ref(error))
4-element Vector{Int64}:
 1
 0
 0
 1
```
"""
function pauli_commutation(a::Int, b::Int)
    return a == 0 || b == 0 || a == b ? 0 : 1
end
function pauli_commutation(a::AbstractVector{Int}, b::AbstractVector{Int})
    return mod(sum(pauli_commutation.(a, b)), 2)
end

"""
    pauli_pow(a::Int, power::Int) -> Int

Return the Pauli raised to the given power.

# Examples
```jldoctest
julia> x_squared = pauli_pow(1, 2)
0

julia> z_inverse = pauli_pow(3, -1)
3
```
"""
function pauli_pow(a::Int, power::Int)
    return mod(power, 2) == 0 ? 0 : a
end

"""
    pauli_product(a::Int, b::Int) -> Int

Return the product of two Paulis.

# Examples
```jldoctest
julia> pauli_product(1, 2) # X.Y -> Z
3

julia> pauli_product.([1, 0, 3, 2], [1, 1, 1, 3]) # (XIZY).(XXXZ) -> (IXYX)
4-element Vector{Int64}:
 0
 1
 2
 1
```
"""
function pauli_product(a::Int, b::Int)
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
    pauli_product_pow(operators, powers)

Return the product of the Pauli operators each raised to the corresponding power, where
`operators` is an iterable of `AbstractVector{Int}` and `powers` is an iterable of `Int`.

Note: The product is evaluated to the length of the shorter of the two iterables.

See also [`pauli_product`](@ref), [`pauli_pow`](@ref).

# Examples
```jldoctest
julia> ops = [[3, 3, 0], [0, 1, 1], [2, 0, 2]];  # ZZI, IXX, YIY

julia> pauli_product_pow(ops, [1, 0, 1])
3-element Vector{Int64}:
 1
 3
 2
```
"""
function pauli_product_pow(operators, powers)
    return reduce(
        (a, b)->pauli_product.(a, b),
        (pauli_pow.(operator, power) for (operator, power) in zip(operators, powers))
    )
end

"""
    pauli_rep_change(pauli::Int) -> Char
    pauli_rep_change(pauli::Char) -> Int

Convert Pauli between alphabetical and numerical representation.

# Examples
```jldoctest
julia> pauli_rep_change(1)
'X': ASCII/Unicode U+0058 (category Lu: Letter, uppercase)

julia> pauli_rep_change('X')
1

julia> pauli_rep_change.((0, 1, 2, 3))
('I', 'X', 'Y', 'Z')

julia> pauli_rep_change.(('I', 'X', 'Y', 'Z'))
(0, 1, 2, 3)
```
"""
function pauli_rep_change(pauli::Int)
    return "IXYZ"[pauli+1]
end
const _pauli_rep_map = Dict(zip("IXYZ", 0:3))
function pauli_rep_change(pauli::Char)
    return _pauli_rep_map[pauli]
end

"""
    pauli_weight(operator::AbstractVector{Int}) -> Int

Return the weight of the Pauli operator, i.e. the number of non-identity elements.

# Examples
```jldoctest
julia> pauli_weight([3, 0, 2, 0, 1])  # ZIYIX
3
```
"""
function pauli_weight(operator::AbstractVector{Int})
    return count(!=(0), operator)
end

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
    pauli_are_independent(operators) -> Bool

Return true if the Pauli operators are independent, or false otherwise, where
`operators` is an iterable of `AbstractVector{Int}`.  Here independence means
that no operator can be expressed as a product of the others.

# Examples
```jldoctest
julia> pauli_are_independent([[3, 3, 0], [0, 3, 3]])  # ZZI, IZZ
true

julia> pauli_are_independent([[3, 3, 0], [0, 3, 3], [3, 0, 3]])  # ZZI, IZZ, ZIZ
false
```
"""
function pauli_are_independent(operators)
    operator_list = collect(operators)  # collect iterable to use length and index access

    num_operators = length(operator_list)
    if num_operators == 0  # empty operator_list is independent
        return true
    end
    num_qubits = length(operator_list[1])
    remaining = collect(1:num_operators)

    for qubit in 1:num_qubits
        paulis = [1,2,3]
        indices = Int[]
        for α in remaining
            if operator_list[α][qubit] in paulis
                setdiff!(paulis, operator_list[α][qubit])
                push!(indices, α)
            end
            if length(paulis) == 1
                break
            end
        end

        remaining = setdiff(remaining, indices)
        for α in remaining
            if operator_list[α][qubit] == 0
                continue
            end
            new_operator = pauli_product.(operator_list[α], operator_list[indices[1]])
            if new_operator[qubit] == 0
                operator_list[α] = deepcopy(new_operator)
                continue
            end
            if length(paulis) == 2
                continue
            end
            new_operator = pauli_product.(operator_list[α], operator_list[indices[2]])
            if new_operator[qubit] == 0
                operator_list[α] = deepcopy(new_operator)
                continue
            end
            new_operator = pauli_product.(operator_list[α], operator_list[indices[1]])
            new_operator = pauli_product.(new_operator, operator_list[indices[2]])
            if new_operator[qubit] == 0
                operator_list[α] = deepcopy(new_operator)
                continue
            end
        end

        for α in remaining
            if operator_list[α] == zeros(Int64, num_qubits)
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
    return mod(sum(pauli_commutation(p1, p2) for (p1, p2) in zip(a, b); init=0), 2)
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

Return the product of two Paulis.  Note that we are ignoring overall minus
signs.

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
    if a == b  # pauli squared is identity
        return 0
    elseif a == 0 || b == 0  # pauli times identity is pauli
        return a + b
    else  # product of distinct paulis is the other pauli
        # X.Y=Z -> 6-(1+2)=3,  X.Z=Y -> 6-(1+3)=2,  Y.Z=X -> 6-(2+3)=1
        return 6 - (a + b)
    end
end

"""
    pauli_product(operators) -> Vector{Int}

Return the product of the Pauli operators, where `operators` is an iterable of
`AbstractVector{Int}`.

See also [`pauli_product_pow`](@ref).

# Examples
```jldoctest
julia> ops = [[3, 3, 0], [0, 1, 1], [2, 0, 2]];  # ZZI, IXX, YIY

julia> pauli_product(ops)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function pauli_product(operators)
    return reduce((a, b)->pauli_product.(a, b), operators)
end

"""
    pauli_product_pow(operators, powers) -> Vector{Int}

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
    pauli_random_operator(n::Int, p::Real, rng::AbstractRNG=GLOBAL_RNG)
        -> AbstractVector{Int}

Return a random Pauli operator on `n` qubits with Paulis applied according to i.i.d.
depolarizing noise with probability `p`.

# Examples
```jldoctest
julia> using Random:MersenneTwister  # use RNG for reproducible example

julia> pauli_random_operator(5, 0.2, MersenneTwister(13))
5-element Vector{Int64}:
 3
 3
 0
 2
 0
```
"""
function pauli_random_operator(n::Int, p::Real, rng::AbstractRNG=GLOBAL_RNG)
    weights = ProbabilityWeights([1 - p, p / 3, p / 3, p / 3])
    return sample(rng, [0, 1, 2, 3], weights, n)
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

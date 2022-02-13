"""
    five_qubit_code() -> SimpleCode

Return the five-qubit code.
"""
function five_qubit_code()
    stabilizers = [[1,3,3,1,0], [0,1,3,3,1], [1,0,1,3,3], [3,1,0,1,3]]
    logicals = [[1,1,1,1,1], [3,3,3,3,3]]
    pure_errors = [[0,1,0,0,0], [0,0,0,0,3], [0,0,3,0,0], [1,0,0,0,0]]
    return SimpleCode("Five qubit code", stabilizers, logicals, pure_errors)
end

"""
    five_qubit_surface_code() -> SimpleCode

Return the five-qubit surface code.
"""
function five_qubit_surface_code()
    stabilizers = [[1,1,1,0,0], [0,0,1,1,1], [3,0,3,3,0], [0,3,3,0,3]]
    logicals = [[1,0,0,1,0], [3,3,0,0,0]]
    pure_errors = [[3,0,0,0,0], [0,0,0,3,0], [1,0,0,0,0], [0,1,0,0,0]]
    return SimpleCode("Five qubit surface code", stabilizers, logicals, pure_errors)
end

"""
    steane_code() -> SimpleCode

Return the seven-qubit Steane code.
"""
function steane_code()
    stabilizers = [[1,0,0,1,0,1,1], [0,1,0,1,1,0,1], [0,0,1,0,1,1,1], [3,0,0,3,0,3,3],
        [0,3,0,3,3,0,3], [0,0,3,0,3,3,3]]
    logicals = [[1,1,1,1,1,1,1], [3,3,3,3,3,3,3]]
    pure_errors = [[3,0,0,0,0,0,0], [0,3,0,0,0,0,0], [0,0,3,0,0,0,0], [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0], [0,0,1,0,0,0,0]]
    return SimpleCode("Steane code", stabilizers, logicals, pure_errors)
end

"""
    random_code(n::Int, k::Int) -> SimpleCode

Return a random code with `n` physical qubits and `k` logicals qubits.

An `ErrorException` is raised if `k` > `n`.
"""
function random_code(n::Int, k::Int)
    0 <= k <= n || error("invalid number of logical qubits!")

    r = n - k

    stabilizer_state = random_stabilizer_state(n)
    stabilizers = stabilizer_state.stabilizers[1:r]
    pure_errors = stabilizer_state.pure_errors[1:r]

    logicals = Vector{Int}[]
    for α in r + 1:n
        push!(logicals, stabilizer_state.stabilizers[α])
        push!(logicals, stabilizer_state.pure_errors[α])
    end

    return SimpleCode("[$n,$k] random code", stabilizers, logicals, pure_errors)
end

"""
    random_stabilizer_state(n::Int) -> SimpleCode

Return a random code with `n` physical qubits and no logicals.
"""
function random_stabilizer_state(n::Int)
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

    logicals = []

    return SimpleCode("$n-qubit random stabilizer state", stabilizers, logicals)
end

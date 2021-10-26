"""
    five_qubit_code()

Returns the five-qubit code as a `SimpleCode`.
"""
function five_qubit_code()

    stabilizers = [[1,3,3,1,0],[0,1,3,3,1],[1,0,1,3,3],[3,1,0,1,3]]
    logicals = [[1,1,1,1,1],[3,3,3,3,3]]
    pure_errors = [[0,1,0,0,0],[0,0,0,0,3],[0,0,3,0,0],[1,0,0,0,0]]

    return SimpleCode(
        "Five qubit code",
        stabilizers,
        logicals,
        pure_errors)

end





"""
    five_qubit_surface_code()

Returns the five-qubit surface code as a `SimpleCode`.
"""
function five_qubit_surface_code()

    stabilizers = [[1,1,1,0,0],[0,0,1,1,1],[3,0,3,3,0],[0,3,3,0,3]]
    logicals = [[1,0,0,1,0],[3,3,0,0,0]]
    pure_errors = [[3, 0, 0, 0, 0],[0, 0, 0, 3, 0],[1, 0, 0, 0, 0],[0, 1, 0, 0, 0]]

    return SimpleCode(
        "Five qubit surface code",
        stabilizers,
        logicals,
        pure_errors)

end




"""
    steane_code()

Returns the seven-qubit Steane code as a `SimpleCode`.
"""
function steane_code()

    stabilizers = [[1,0,0,1,0,1,1],[0,1,0,1,1,0,1],[0,0,1,0,1,1,1],
    [3,0,0,3,0,3,3],[0,3,0,3,3,0,3],[0,0,3,0,3,3,3]]
    logicals = [[1,1,1,1,1,1,1],[3,3,3,3,3,3,3]]
    pure_errors = [[3,0,0,0,0,0,0],[0,3,0,0,0,0,0],[0,0,3,0,0,0,0],
    [1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0]]

    return SimpleCode(
        "Steane code",
        stabilizers,
        logicals,
        pure_errors)

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

    logicals = []

    return SimpleCode(
        "random stabilizer state",
        stabilizers,
        logicals)
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

"""
    min_weight_brute_force(code,syndrome,
    error_prob,error_model)

Returns an error which is consistent with the `syndrome` for this code
and which has minimum weight.  No effort has been made to make this fast.
(Though the problem is exponential time anyway!)
"""
function min_weight_brute_force(
        code::QuantumCode,
        syndrome::Array{Int64,1};
        error_model="depolarizing")

    g = code.stabilizers
    r = length(g)
    n = num_qubits(code)

    if r != length(syndrome)
        error("number of syndrome bits does not match number of stabilizers!")
    end

    error_guess = zeros(Int64,n)
    error_weight = n + 1
    for α in 0:4^n-1
        operator = digits!(zeros(Int64,n),α,base = 4)
        w = pauli_weight(operator)
        if (find_syndrome(code,operator) == syndrome &&
                w <= error_weight)
            error_weight = w
            error_guess = operator
        end
    end

    return error_guess
end





"""
    do_nothing_decoder(code,syndrome,error_prob,error_model)

Returns a correction which does nothing!
"""
function do_nothing_decoder(code::QuantumCode,syndrome::Array{Int64,1},
        error_prob)

    n = num_qubits(code)

    return zeros(Int64,n)
end





"""
    monte_carlo_simulation(code,probabilities,N,decoder)

Randomly generates errors, decodes and then evaluates success
probabilities for a `code` given a range of physical error
`probabilities` each for `N` Monte Carlo samples using `decoder`.
"""
function monte_carlo_simulation(
        code::QuantumCode,
        probabilities::Array{Float64,1},
        N::Int64,
        decoder)

    successes = Float64[]

    for p in probabilities
        success = 0
        for seed in 1:N
            initial_error = random_pauli_error(num_qubits(code),p,seed)
            syndrome = find_syndrome(code,initial_error)

            correction = decoder(code,syndrome,p)
            effect_on_code = pauli_product.(initial_error,correction)

            if (pauli_commutation.(Ref(effect_on_code),code.logicals)
                == zeros(Int64,length(code.logicals)))
                success += 1
            end
        end

        push!(successes,success)
    end

    return successes./N
end

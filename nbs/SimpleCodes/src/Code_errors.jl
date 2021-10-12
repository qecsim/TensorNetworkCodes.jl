import StatsBase
using Random




"""
    random_pauli_error(n::Int64,error_prob::Float64,seed::Int64)

Returns random error according to i.i.d. depolarizing noise on `n` qubits.
"""
function random_pauli_error(n::Int64,error_prob::Float64,seed::Int64)
    rng = MersenneTwister(seed)
    values = [0, 1, 2, 3]

    probabilities = [1-error_prob, error_prob/3, error_prob/3, error_prob/3]
    w = StatsBase.Weights(probabilities)
    return [StatsBase.sample(rng,values,w) for α in 1:n]
end



function random_pauli_error(n::Int64,error_prob::Float64)
    rng = RandomDevice()
    values = [0, 1, 2, 3]

    probabilities = [1-error_prob, error_prob/3, error_prob/3, error_prob/3]
    w = StatsBase.Weights(probabilities)
    return [StatsBase.sample(rng,values,w) for α in 1:n]
end
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
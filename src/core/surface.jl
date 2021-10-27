"""
    _diamond_lattice_code(code_matrix::Matrix{TensorNetworkCode})

Takes a matrix of `TensorNetworkCode`s and returns a big `TensorNetworkCode` with all
constituent codes contracted in a diamond lattice.
"""
function _diamond_lattice_code(code_matrix::Matrix{TensorNetworkCode})

    L,W = size(code_matrix)

    # Coordinate shifts for each code in code_matrix
    coord_shifts =
    [(isodd(l) ? [0,0] : [6,0]) + [12*(w-1),6*(l-1)]
        for l in 1:L, w in 1:W]

    # Starting coords for each node (first is virtual & the rest are
    # physical qubits) in a seed code
    coords = [[0,0],[3,-3],[3,3],[0,0.5],[-3,-3],[-3,3]]

    output_code = code_matrix[1,1]
    set_coords!(output_code,coords)

    # contract on the rest of the codes in code_matrix
    for α in 1:L,β in 1:W
        if (α == 1 && β == 1) || (iseven(α) && β == W)
            continue
        end

        next_seed_code = code_matrix[α,β]
        set_coords!(next_seed_code,coords)
        shift = coord_shifts[α,β]
        shift_coords!(next_seed_code,shift)

        output_code = contract_by_coords(output_code,next_seed_code)

        if num_qubits(output_code) == 0
            return output_code
        end
    end

    return output_code
end



_diamond_lattice_code(code_matrix::Matrix{SimpleCode}) =
_diamond_lattice_code(TensorNetworkCode.(code_matrix))





"""
    surface_code_bulk(L::Int64)

Generates the bulk `TensorNetworkCode` for the surface code.
"""
function _surface_code_bulk(L::Int64)

    # Surface code with 5 qubits
    small_surface = five_qubit_surface_code()

    # Other constituent codes are made by choosing certain logical
    # operators to be stabilizers, resulting in stabilizers states
    x_small_surface = gauge(small_surface, 1, 1)
    z_small_surface = gauge(small_surface, 1, 3)

    # this one is rotated:
    z_small_surface_rotated = permute(z_small_surface,[2,5,3,1,4])

    code_matrix = [isodd(i) ? x_small_surface : z_small_surface_rotated
        for i in 1:2*L-1, j in 1:L]

    code_matrix[1,1] = small_surface # this is the only one with a logical
    for i in 2:2*L-1
        if isodd(i)
            code_matrix[i,1] = z_small_surface
        end
    end

    return _diamond_lattice_code(code_matrix)
end





"""
    boundary_fusion(
        L::Int64,
        bulk_code_function=surface_code_bulk)

To create the surface code, we need to get the right boundary
conditions, which this code does.
"""
function _boundary_fusion(L::Int64,bulk_code::TensorNetworkCode)


    # Repitition codes for boundary contraction
    x_repitition_zero = SimpleCode(
        "X repitition code",
        [[3,3,0],[0,3,3],[1,1,1]],
        [],
        [[1,0,0],[0,0,1],[3,0,0]])
    z_repitition_zero = SimpleCode(
        "Z repitition code",
        [[1,1,0],[0,1,1],[3,3,3]],
        [],
        [[3,0,0],[0,0,3],[1,0,0]])
    z_rep = TensorNetworkCode(z_repitition_zero)
    x_rep = TensorNetworkCode(x_repitition_zero)


    # depending on the side (top, left, etc), we choose different coordinates
    left_z_rep = deepcopy(z_rep)
    coords = [[0,0],[-3,0],[-3,-3],[-3,3]]
    set_coords!(left_z_rep,coords)

    right_z_rep = z_rep
    coords = [[0,0],[3,0],[3,-3],[3,3]]
    set_coords!(right_z_rep,coords)

    top_x_rep = deepcopy(x_rep)
    coords = [[0,0],[0,-3],[-3,-3],[3,-3]]
    set_coords!(top_x_rep,coords)

    bottom_x_rep = x_rep
    coords = [[0,0],[0,3],[-3,3],[3,3]]
    set_coords!(bottom_x_rep,coords)


    # Now we can put the code together
    output = bulk_code
    if num_qubits(output) == 0
        return output
    end


    # Let's start fusing boundary codes, starting with the top row
    code = top_x_rep
    for l in 1:L-1
        l==1 ? shift = [6,0] : shift = [12,0]
        shift_coords!(code,shift)
        output = contract_by_coords(output,code)
    end

    code = bottom_x_rep
    for l in 1:L-1
        l == 1 ? shift = [6,12*(L-1)] : shift = [12,0]
        shift_coords!(code, shift)
        output = contract_by_coords(output,code)
    end

    code = left_z_rep
    for l in 1:L-1
        l==1 ? shift = [0,6] : shift = [0,12]
        shift_coords!(code, shift)
        output = contract_by_coords(output,code)
    end

    code = right_z_rep
    for l in 1:L-1
        l==1 ? shift = [12*(L-1),6] : shift = [0,12]
        shift_coords!(code, shift)
        output = contract_by_coords(output,code)
    end

    return output
end





"""
    surface_code(L::Int64)

Returns an (L+1) x L surface code `TensorNetworkCode`.

# Examples
```jldoctest
julia> surf = surface_code(2);

julia> dist = find_distance_logicals(surf)[1] # code distance
3
```
"""
function surface_code(L::Int64)

    bulk_code = _surface_code_bulk(L)
    return _boundary_fusion(L,bulk_code)
end



"""
    almost_rotated_surface_code(
    L::Int64,
    input_seed_code::SimpleCode,
    input_coords::Array{Int64}) -> TensorNetworkCode

Generates a `TensorNetworkCode` that is almost the rotated surface
code.  Requires the `input_seed_code` to have five qubits and to be
in the bulk (`input_coords` with components not equal to 1 or L).

# Examples
```jldoctest
julia> code = almost_rotated_surface_code(3,five_qubit_code(),[2,2]);

julia> verify_code(code) # try this out and check we get a real code
true
```
"""
function almost_rotated_surface_code(
    L::Int64,
    input_seed_code::SimpleCode,
    input_coords::Array{Int64})

"""
The codes and labels below refer to seed codes laid out in the
following pattern.

3x3:
A  B  C
D  S  F
G  H  I

5x5:
A  B  Br B  C
D  E  Er E  F
Dr Er S  Er Fr
D  E  Er E  F
G  H  Hr H  I

S is a five-qubit surface code with a single logical.
E is a five qubit code with a logical as a stabilizer.
Er is the same as E with Hadamards applied to every qubit.
The rest are detailed below.
"""

if iseven(L)
    error("This only works for odd L!")
end

code_dictionary = Dict{String,TensorNetworkCode}()

# A
stabilizers = [[1,1,1],[0,3,3],[3,0,3]]
logicals = []
code = SimpleCode("A",stabilizers,logicals)
A = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[1,0],[0.3,0.3]]
set_coords!(A,coordinates)
code_dictionary["A"] = A

# B
stabilizers = [[1,1,0,1],[3,0,0,3],[3,0,3,3],[0,3,0,3]]
logicals = []
code = SimpleCode("B",stabilizers,logicals)
B = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[1,0],[0.3,0.3]]
set_coords!(B,coordinates)
code_dictionary["B"] = B

# Br
stabilizers = [[3,3,0,3],[0,1,1,1],[1,0,1,1],[0,0,3,3]]
logicals = []
code = SimpleCode("Br",stabilizers,logicals)
Br = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[1,0],[0.3,0.3]]
set_coords!(Br,coordinates)
code_dictionary["Br"] = Br

# C
stabilizers = [[3,3,3],[0,1,1],[1,0,1]]
logicals = []
code = SimpleCode("C",stabilizers,logicals)
C = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[0.3,0.3]]
set_coords!(C,coordinates)
code_dictionary["C"] = C

# D
stabilizers = [[0,0,1,1],[0,1,1,1],[3,0,3,3],[1,0,0,1]]
logicals = []
code = SimpleCode("D",stabilizers,logicals)
D = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[0,-1],[1,0],[0.3,0.3]]
set_coords!(D,coordinates)
code_dictionary["D"] = D

# Dr
stabilizers = [[3,3,0,3],[0,3,3,3],[1,0,1,1],[0,1,0,1]]
logicals = []
code = SimpleCode("Dr",stabilizers,logicals)
Dr = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[0,-1],[1,0],[0.3,0.3]]
set_coords!(Dr,coordinates)
code_dictionary["Dr"] = Dr

# E
stabilizers = [[0,1,1,0,1],[1,0,0,1,1],[3,3,0,0,3],[0,0,3,3,3],[1,1,0,0,0]]
logicals = []
code = SimpleCode("E",stabilizers,logicals)
E = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(E,coordinates)
code_dictionary["E"] = E

# Er
stabilizers = [[0,3,3,0,3],[3,0,0,3,3],[1,1,0,0,1],[0,0,1,1,1],[3,3,0,0,0]]
logicals = []
code = SimpleCode("Er",stabilizers,logicals)
Er = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(Er,coordinates)
code_dictionary["Er"] = Er

# F
stabilizers = [[0,3,3,3],[1,1,0,1],[0,1,0,1],[0,0,1,1]]
logicals = []
code = SimpleCode("F",stabilizers,logicals)
F = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[0,-1],[0.3,0.3]]
set_coords!(F,coordinates)
code_dictionary["F"] = F

# Fr
stabilizers = [[0,1,1,1],[3,3,0,3],[3,0,3,3],[1,0,0,1]]
logicals = []
code = SimpleCode("Fr",stabilizers,logicals)
Fr = TensorNetworkCode(code)
coordinates = [[0,0],[0,1],[-1,0],[0,-1],[0.3,0.3]]
set_coords!(Fr,coordinates)
code_dictionary["Fr"] = Fr

# G
stabilizers = [[3,3,3],[1,0,1],[0,1,1]]
logicals = []
code = SimpleCode("G",stabilizers,logicals)
G = TensorNetworkCode(code)
coordinates = [[0,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(G,coordinates)
code_dictionary["G"] = G

# H
stabilizers = [[3,3,0,3],[0,1,1,1],[0,3,0,3],[0,0,3,3]]
logicals = []
code = SimpleCode("H",stabilizers,logicals)
H = TensorNetworkCode(code)
coordinates = [[0,0],[-1,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(H,coordinates)
code_dictionary["H"] = H

# Hr
stabilizers = [[1,1,0,1],[0,3,3,3],[1,0,1,1],[3,0,0,3]]
logicals = []
code = SimpleCode("Hr",stabilizers,logicals)
Hr = TensorNetworkCode(code)
coordinates = [[0,0],[-1,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(Hr,coordinates)
code_dictionary["Hr"] = Hr

# I
stabilizers = [[1,1,1],[3,0,3],[0,3,3]]
logicals = []
code = SimpleCode("I",stabilizers,logicals)
I = TensorNetworkCode(code)
coordinates = [[0,0],[-1,0],[0,-1],[0.3,0.3]]
set_coords!(I,coordinates)
code_dictionary["I"] = I


# Assign code tensor labels to a matrix
code_label_matrix = fill("", L, L)
for i in 1:L, j in 1:L
    if iseven(i+j)
        code_label_matrix[i,j] = "E"
    else
        code_label_matrix[i,j] = "Er"
    end
end

for j in 1:L
    if iseven(j)
        code_label_matrix[1,j] = "B"
    else
        code_label_matrix[1,j] = "Br"
    end
end

for j in 1:L
    if iseven(j)
        code_label_matrix[L,j] = "H"
    else
        code_label_matrix[L,j] = "Hr"
    end
end

for i in 1:L
    if iseven(i)
        code_label_matrix[i,1] = "D"
    else
        code_label_matrix[i,1] = "Dr"
    end
end

for i in 1:L
    if iseven(i)
        code_label_matrix[i,L] = "F"
    else
        code_label_matrix[i,L] = "Fr"
    end
end

# Take care of the corners
code_label_matrix[1,1] = "A"
code_label_matrix[1,L] = "C"
code_label_matrix[L,1] = "G"
code_label_matrix[L,L] = "I"


# Include input code
if num_qubits(input_seed_code) != 5
    error("input_seed_code does not have five physical qubits!")
elseif 1 in input_coords || L in input_coords
    error("location of input_seed_code needs to be in the bulk!")
end

input_tn_code = TensorNetworkCode(input_seed_code)
coordinates = [[0,0],[0,1],[-1,0],[0,-1],[1,0],[0.3,0.3]]
set_coords!(input_tn_code,coordinates)
code_dictionary["input"] = input_tn_code

# Add to the code label matrix
code_label_matrix[input_coords...] = "input"


output_code = deepcopy(code_dictionary[code_label_matrix[1,1]])
for i in 1:L
    for j in 1:L
        if i==j && i==1
            continue
        end
        code = deepcopy(code_dictionary[code_label_matrix[i,j]])
        shift = [Float64(j-1)*2,Float64(i-1)*2]
        shift_coords!(code,shift)
        output_code = contract_by_coords(output_code,code)
    end
end

return output_code
end





"""
    rotated_surface_code(L::Int64) -> TensorNetworkCode

Generates the LxL rotated surface code as a `TensorNetworkCode`.

# Examples
```jldoctest
julia> rot = rotated_surface_code(3);

julia> pauli_weight(rot.stabilizers[2]) # weight of a plaquette stabilizer
4
```
"""
function rotated_surface_code(L::Int64)

stabilizers = [[0,1,1,0,1],[1,0,0,1,1],[3,3,0,0,3],[0,0,3,3,3]]
logicals = [[3,0,3,0,3],[0,1,0,1,1]]
code = SimpleCode("central",stabilizers,logicals)

# Add to the code label matrix
centre = ceil(Int,L/2)
input_coords = [centre,centre]

return almost_rotated_surface_code(L,code,input_coords)
end

"""
    diamond_lattice_code(codes::Matrix{TNCode})

Takes a matrix of `TNCodes` and returns a big `TNCode` with all
constituent codes contracted in a diamond lattice.
"""
function diamond_lattice_code(codes::Matrix{TNCode})

    L,W = size(codes)
    coords = [[0,0],[3,-3],[3,3],[0,0.5],[-3,-3],[-3,3]]

    output = deepcopy(codes[1,1])
    set_coords!(output,coords)


    odd_first(l) = [0,+6*l]
    even_first(l) = [6,6+6*(l-1)]
    vec = [12,0]
    sws = [isodd(l) ? odd_first(l-1) .+ vec*(w-1) :
    even_first(l-1) .+ vec*(w-1) for l in 1:L, w in 1:W]


    for α in 1:L,β in 1:W
        if (α == 1 && β == 1) || (iseven(α) && β == W)
            continue
        end

        next_code = codes[α,β]
        set_coords!(next_code,coords)
        shift = sws[α,β]
        shift_coords!(next_code,shift)

        output = combine_by_coordinates(output,next_code)

        if size(output) == 0
            return output
        end
    end

    return output
end



diamond_lattice_code(codes::Matrix{SimpleCode}) =
diamond_lattice_code(TNCode.(codes))





"""
    surface_code_bulk(L::Int64)

Generates the bulk `TNCode` for the surface code.
"""
function surface_code_bulk(L::Int64)


    # Surface code with 5 qubits
    small_surface = SimpleCode("Small surface code",
    [[1,1,1,0,0],[0,0,1,1,1],[3,0,3,3,0],[0,3,3,0,3]],
    [[1,0,0,1,0],[3,3,0,0,0]],
    [[3, 0, 0, 0, 0],[0, 0, 0, 3, 0],[1, 0, 0, 0, 0],[0, 1, 0, 0, 0]])

    # Other constituent codes are made by choosing certain logical
    # operators to be stabilizers, resulting in stabilizers states
    x_small_surface = gauge_code(small_surface,[[1,0]],[1])
    z_small_surface = gauge_code(small_surface,[[0,1]],[1])

    z_small_surface_rotated = TensorNetworkCodes.permute(z_small_surface,[2,5,3,1,4])  # this is rotated
    # z_small_surface_rotated = SimpleCode("Rotated small surface code",
    # z_small_surface_rot.stabilizers,
    # z_small_surface_rot.logicals,
    # z_small_surface_rot.pure_errors)



    codes = [isodd(i) ? x_small_surface : z_small_surface_rotated
        for i in 1:2*L-1, j in 1:L]

    codes[1,1] = small_surface # this is the only one with a logical
    for i in 2:2*L-1
        if isodd(i)
            codes[i,1] = z_small_surface
        end
    end


    return diamond_lattice_code(codes)
end





"""

"""
function almost_surface_code_bulk(
        L::Int64,
        input_code::SimpleCode,
        location::Array{Int64,1}=[1,1])


    # Surface code with 5 qubits
    small_surface = SimpleCode("Small surface code",
    [[1,1,1,0,0],[0,0,1,1,1],[3,0,3,3,0],[0,3,3,0,3]],
    [[1,0,0,1,0],[3,3,0,0,0]],
    [[3, 0, 0, 0, 0],[0, 0, 0, 3, 0],[1, 0, 0, 0, 0],[0, 1, 0, 0, 0]])

    # Other constituent codes are made by choosing certain logical
    # operators to be stabilizers, resulting in stabilizers states
    x_small_surface = gauge_code(small_surface,[[1,0]],[1])
    z_small_surface = gauge_code(small_surface,[[0,1]],[1])
    z_small_surface_rotated = TensorNetworkCodes.permute(z_small_surface,[2,5,3,1,4])  # this is rotated


    codes = [isodd(i) ? x_small_surface : z_small_surface_rotated
        for i in 1:2*L-1, j in 1:L]

    codes[1,1] = small_surface # this is the only one with a logical
    for i in 2:2*L-1
        if isodd(i)
            codes[i,1] = z_small_surface
        end
    end


    codes[location...] = input_code

    return diamond_lattice_code(codes)
end





"""
    boundary_fusion(
        L::Int64,
        bulk_code_function=surface_code_bulk)

To create the surface code, we need to get the right boundary
conditions, which this code does.
"""
function boundary_fusion(
        L::Int64,
        bulk_code::TNCode)


    # Repitition codes
    x_repitition_zero = SimpleCode("X repitition code",[[3,3,0],[0,3,3],[1,1,1]],[],
    [[1,0,0],[0,0,1],[3,0,0]])
    z_repitition_zero = SimpleCode("Z repitition code",[[1,1,0],[0,1,1],[3,3,3]],[],
    [[3,0,0],[0,0,3],[1,0,0]])
    z_rep = TNCode(z_repitition_zero)
    x_rep = TNCode(x_repitition_zero)


    # depending on the side (top,left, etc), we choose different coordinates
    left_z_rep = deepcopy(z_rep)
    coords = [[0,0],[-3,0],[-3,-3],[-3,3]]
    set_coords!(left_z_rep,coords)

    right_z_rep = deepcopy(z_rep)
    coords = [[0,0],[3,0],[3,-3],[3,3]]
    set_coords!(right_z_rep,coords)

    top_x_rep = deepcopy(x_rep)
    coords = [[0,0],[0,-3],[-3,-3],[3,-3]]
    set_coords!(top_x_rep,coords)

    bottom_x_rep = deepcopy(x_rep)
    coords = [[0,0],[0,3],[-3,3],[3,3]]
    set_coords!(bottom_x_rep,coords)


    # Now we can put the code together
    output = bulk_code
    if size(output) == 0
        return output
    end


    # Let's start fusing boundary codes
    code = SimpleCode()
    for l in 1:L-1
        if l==1
            code = deepcopy(top_x_rep)
            shift_coords!(code, [6,0])
            output = combine_by_coordinates(output,code)
        else
            shift_coords!(code, [12,0])
            output = combine_by_coordinates(output,code)
        end
    end


    for l in 1:L-1
        if l==1
            code = deepcopy(bottom_x_rep)
            shift_coords!(code, [6,12*(L-1)])
            output = combine_by_coordinates(output,code)
        else
            shift_coords!(code, [12,0])
            output = combine_by_coordinates(output,code)
        end
    end


    for l in 1:L-1
        if l==1
            code = deepcopy(left_z_rep)
            shift_coords!(code, [0,6])
            output = combine_by_coordinates(output,code)
        else
            shift_coords!(code, [0,12])
            output = combine_by_coordinates(output,code)
        end
    end


    for l in 1:L-1
        if l==1
            code = deepcopy(right_z_rep)
            shift_coords!(code, [12*(L-1),6])
            output = combine_by_coordinates(output,code)
        else
            shift_coords!(code, [0,12])
            output = combine_by_coordinates(output,code)
        end
    end


    return output
end





"""
    surface_code(L::Int64)

Returns an LxL surface code `TNCode`.
"""
function surface_code(L::Int64)

    bulk_code = surface_code_bulk(L)
    return boundary_fusion(L,bulk_code)
end




"""
    almost_surface_code(L,input_code,location=[1,1])

Almost a surface code but with a custom `input_code` at `location`.
"""
function almost_surface_code(
        L::Int64,
        input_code::SimpleCode,
        location::Array{Int64,1}=[1,1])

    bulk_code = almost_surface_code_bulk(L,input_code,location)

    return boundary_fusion(L,bulk_code)
end





"""
    fully_random_code(L)

Returns a `TNCode` with the same tensor network as the surface code
but with the bulk tensors replaced by random code tensors.
"""
function fully_random_code(L::Int64)

    codes = [random_stabilizer_state(5) for i in 1:2*L-1, j in 1:L]
    centre = ceil(Int,L/2)
    codes[centre,centre] = random_code(5,1) # this is the only one with a logical

    bulk_code =  diamond_lattice_code(codes)

    return boundary_fusion(L,bulk_code)
end







function checkerboard_code(
        L::Int64,
        x::TNCode,
        z::TNCode,
        centre_code::TNCode)

    codes = [isodd(i) ? z : x
        for i in 1:2*L-1, j in 1:L]

    centre = ceil(Int,L/2)
    codes[centre,centre] = centre_code # this is the only one with a logical

    bulk_code = diamond_lattice_code(codes)

    return boundary_fusion(L,bulk_code)
end





function funny_code(L::Int64)

    # The big code is built from three smaller codes (all based on the usual
    # five-qubit code).
    stabilizers = [[0,0,2,2,1],[1,0,3,0,3],[0,2,1,0,2],[2,0,0,1,2],[3,3,0,0,1]]
    pure_errors = generate_pure_errors(stabilizers)
    logicals = []
    z = TNCode(SimpleCode("z",stabilizers,logicals,pure_errors))

    stabilizers = [[0,0,2,2,3],[2,0,0,3,2],[1,1,0,0,3],[0,2,3,0,2],[0,3,0,1,1]]
    pure_errors = generate_pure_errors(stabilizers)
    x = TNCode(SimpleCode("x",stabilizers,logicals,pure_errors))


    # Finally the only code with a logical qubit is `centre` below.
    stabilizers = [[0,0,2,2,1],[0,2,1,0,2],[2,0,0,1,2],[3,3,0,0,1]]
    pure_errors = generate_pure_errors(stabilizers)
    logicals = [[1,0,3,0,3],[3,0,0,2,0]]
    centre_code = TNCode(SimpleCode("centre",stabilizers,logicals,pure_errors))


    return checkerboard_code(L,x,z,centre_code)
end

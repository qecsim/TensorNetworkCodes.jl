"""
    code_to_tensor(code)

Returns a tensor (array) describing the logical cosets of the code.
"""
function code_to_tensor(code::QuantumCode)
    
    alt_code = purify_code(code)
    n = size(alt_code)
    g = alt_code.stabilizers
    
    dims = fill(4,n)
    tensor = zeros(Float64,dims...)
    
    for α in 0:2^n-1
        powers = digits!(zeros(Int64,n),α,base = 2)
        operator = SimpleCodes.product_with_powers(g,powers)        
        operator = operator .+ 1      

        index = CartesianIndex(operator...)
        tensor[index] = 1        
    end
    
    return tensor   
end





"""
    code_to_Itensor(code,indices)

Returns an `ITensor` with indices describing the logical cosets of the code.
"""
function code_to_Itensor(code::QuantumCode,
        logical_indices::Array{Index{Int64},1},
        physical_indices::Array{Index{Int64},1})

    n = size(code)
    K = length(code.logicals)
    if K/2 != length(logical_indices) || n != length(physical_indices)
        error("incorrect number of indices!")
    end
    
    tensor = code_to_tensor(code)
    indices = vcat(logical_indices,physical_indices)
    
    return ITensor(tensor,indices...) 
end





"""
    identity_coset(tensor)

Given an `ITensor` descibing a code, return the `ITensor` describing
only the identity coset, i.e., the stabilizer group.
"""
function identity_coset(tensor::ITensor)
    indices = inds(tensor)
    logical_indices = [χ for χ in indices if hastags(χ,"logical")]
    
    for index in logical_indices
        ρ = ITensor(index)
        ρ[index => 1] = 1
        tensor = tensor * ρ
    end
    
    return tensor
end




            
"""
    all_cosets(tensor)

Given an `ITensor` descibing a code, return the `ITensor` describing
all the cosets by summing over all logical indices.
"""
function all_cosets(tensor::ITensor)
    indices = inds(tensor)
    logical_indices = [χ for χ in indices if hastags(χ,"logical")]
    
    for index in logical_indices
        ρ = ITensor(ones(dim(index)),index)
        tensor = tensor * ρ
    end
    
    return tensor
end

       
                        
                        
                        
# Extends the method for `SimpleCodes`.                    
function gauge_code(
        code::TNCode,
        logical_power_list::Array{Array{Int64,1},1},
        which_logicals::Array{Int64,1})
    
    new_code = SimpleCode(code)
    new_code = gauge_code(
        new_code,
        logical_power_list,
        which_logicals)
    
    return TNCode(
        new_code,
        code.code_graph,
        code.seed_codes)
end
"""
    hello(name)

Say hello.
"""
function hello(name)
    return "Hello $(name)"
end

function hello(code::StabilizerCode)
    return "Hello $(label(code)) code"
end

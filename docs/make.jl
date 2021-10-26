using TensorNetworkCodes
using Documenter

DocMeta.setdocmeta!(TensorNetworkCodes, :DocTestSetup, :(using TensorNetworkCodes); recursive=true)

makedocs(;
    modules=[TensorNetworkCodes],
    authors="Terry Farrelly, David K. Tuckett",
    repo="https://github.com/dkt29/TensorNetworkCodes.jl/blob/{commit}{path}#{line}",
    sitename="TensorNetworkCodes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dkt29.github.io/TensorNetworkCodes.jl",
        assets=String[],
    ),
    pages=[
        "Overview" => "index.md",
        "API" => [
            "api/types.md",
            "api/pauli_functions.md",
            "api/code_functions.md",
            "api/code_graph_functions.md",
            "api/contraction_functions.md",
            "api/plotting_functions.md",
            "api/example_codes.md",
            "api/qecsim_adaptors.md",
            "api/index.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/dkt29/TensorNetworkCodes.jl",
    devbranch="main",
)

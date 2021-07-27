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
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dkt29/TensorNetworkCodes.jl",
    devbranch="main",
)

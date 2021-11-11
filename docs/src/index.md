```@meta
CurrentModule = TensorNetworkCodes
```

# TensorNetworkCodes.jl

## Introduction

[TensorNetworkCodes.jl](https://github.com/qecsim/TensorNetworkCodes.jl) is a Julia library
developed to support the following research:

* T. Farrelly, D. K. Tuckett, T. M. Stace, _Local tensor-network codes_,
  [arXiv:2109.11996](https://arxiv.org/abs/2109.11996), (2021).

## Installation

_TensorNetworkCodes.jl_ is installed using the Julia package manager
[Pkg](https://pkgdocs.julialang.org/):

```julia
pkg> # Press ']' to enter the Pkg REPL mode
pkg> add https://github.com/qecsim/Qecsim.jl  # NOT REQUIRED WHEN QECSIM IS REGISTERED
pkg> add https://github.com/qecsim/TensorNetworkCodes.jl

```
or
```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/qecsim/Qecsim.jl")  # NOT REQUIRED WHEN QECSIM IS REGISTERED
julia> Pkg.add(url="https://github.com/qecsim/TensorNetworkCodes.jl")

```

## Demos

_TODO: replace with nbviewer links (works better) when TensorNetworkCodes.jl repo is public_

The following demos correspond to results included in
[arXiv:2109.11996](https://arxiv.org/abs/2109.11996):

* [Small example tensor-network codes and code distances]
  (https://github.com/qecsim/TensorNetworkCodes.jl/blob/main/nbs/Small_examples_with_distance.ipynb)

* [The 19 qubit colour code as a tensor-network code]
  (https://github.com/qecsim/TensorNetworkCodes.jl/blob/main/nbs/Colour_code.ipynb)

* [The modified surface code]
  (https://github.com/qecsim/TensorNetworkCodes.jl/blob/main/nbs/Modified_surface_code_example.ipynb)

## Citing

Please cite _TensorNetworkCodes.jl_ if you use it in your research.

A suitable BibTeX entry is:

    @article{Farrelly_LocalTNCodes_2021,
        title = {Local tensor-network codes},
        author = {Farrelly, Terry and Tuckett, David K. and Stace, Thomas M.},
        year = {2021},
        archiveprefix = {arXiv},
        eprint = {2109.11996},
        url = {https://arxiv.org/abs/2109.11996},
    }

Similarly, please cite [Qecsim.jl](https://github.com/qecsim/Qecsim.jl) if you use its
features in your research, see
[Qecsim.jl Documentation](https://qecsim.github.io/Qecsim.jl/) for details.

## License

_TensorNetworkCodes.jl_ is released under the BSD 3-Clause license, see
[LICENSE](https://github.com/qecsim/TensorNetworkCodes.jl/blob/main/LICENSE).

## Links

* Source code: <https://github.com/qecsim/TensorNetworkCodes.jl>
* Documentation: <https://qecsim.github.io/TensorNetworkCodes.jl>
* Contact: _TODO: add contact email_

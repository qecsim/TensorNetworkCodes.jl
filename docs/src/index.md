```@meta
CurrentModule = TensorNetworkCodes
```

# TensorNetworkCodes.jl

## Introduction

_TensorNetworkCodes.jl_ is a Julia library developed to support the following research:

* T. Farrelly, D. K. Tuckett, T. M. Stace, _Local tensor-network codes_,
  [arXiv:2109.11996](https://arxiv.org/abs/2109.11996), (2021).

## Installation

_TODO: this will be simple when TensorNetworkCodes.jl and Qecsim.jl are publicly released_

## Demos

_TODO: replace with nbviewer links (works better) when TensorNetworkCodes.jl is publicly released_

The following demos correspond to results included in
[arXiv:2109.11996](https://arxiv.org/abs/2109.11996):

* [Small example tensor-network codes and code distances]
  (https://github.com/dkt29/TensorNetworkCodes.jl/blob/main/nbs/Small_examples_with_distance.ipynb)

* [The 19 qubit colour code as a tensor-network code]
  (https://github.com/dkt29/TensorNetworkCodes.jl/blob/main/nbs/Colour_code.ipynb)

* [The modified surface code]
  (https://github.com/dkt29/TensorNetworkCodes.jl/blob/main/nbs/Modified_surface_code_example.ipynb)

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

Similarly, please cite [Qecsim.jl](https://github.com/dkt29/Qecsim.jl) if you use its
features in your research, see
[Qecsim.jl Documentation](https://dkt29.github.io/Qecsim.jl/) for details.

## License

_TensorNetworkCodes.jl_ is released under the BSD 3-Clause license, see
[LICENSE](https://github.com/dkt29/TensorNetworkCodes.jl/blob/main/LICENSE).

### API
```@contents
Pages = [
    "api/types.md",
    "api/pauli_functions.md",
    "api/code_functions.md",
    "api/code_graph_functions.md",
    "api/plotting_functions.md",
    "api/example_codes.md",
    "api/itensors_functions.md",
    "api/TNDecode.md",
    "api/TNDistance.md",
    "api/QecsimAdaptors.md",
    "api/index.md",
]
Depth = 1
```

## Links

* Source code: <https://github.com/dkt29/TensorNetworkCodes.jl>
* Documentation: <https://dkt29.github.io/TensorNetworkCodes.jl/>
* Contact: _TODO: add contact email_

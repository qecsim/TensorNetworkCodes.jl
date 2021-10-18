# Pauli functions

Single-qubit Pauli operators ``I``, ``X``, ``Y``, ``Z`` are represented by
integers 0, 1, 2, 3, respectively. Correspondingly, multi-qubit Pauli operators
are represented by vectors of integers. For example, `[1, 0, 3, 2]` represents
the Pauli operator ``X \otimes I \otimes Z \otimes Y``.

```@docs
pauli_product
do_they_commute
pauli_rep_change
product_with_powers
weight
are_they_independent
```
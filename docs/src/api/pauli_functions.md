# Pauli functions

Pauli functions for evaluating products, commutations, weights and independence.
Single-qubit Pauli operators ``I``, ``X``, ``Y``, ``Z`` are represented by integers
0, 1, 2, 3, respectively. Correspondingly, multi-qubit Pauli operators are represented by
vectors of integers; for example, ``X ⊗ I ⊗ Z ⊗ Y`` is represented by `[1, 0, 3, 2]`.

```@docs
pauli_are_commuting
pauli_are_independent
pauli_commutation
pauli_pow
pauli_product
pauli_product_pow
pauli_rep_change
pauli_weight
```
# MomentHierarchy

This package implements JuMP models for the moment hierarchy in polynomial optimization.

## TO DO

### dual decomposition 

* parallelize build ?

#### subproblems

* distribute each constraint to all supporting variable subsets ? (so far, yes)

### primal-dual decomposition

* implement ?
* if so, move intersecting code with dual decomposition to decomposition_tools.jl

### sparsity 

* Term Sparsity
* tests for Correlative Sparsity

### general

* normalize problem for numerical stability ?
* create Type to dispatch relaxation models and share code ?

### certify

* eigen values computation: https://discourse.julialang.org/t/computation-of-eigenvalues-of-a-sparse-matrix-that-is-given-by-a-file-in-coo-format/75214
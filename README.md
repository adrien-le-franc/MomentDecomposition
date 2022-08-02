# MomentHierarchy

This package implements JuMP models for the moment hierarchy in polynomial optimization.

## TO DO

### dual decomposition 

* subgradient method dependency : choose a reliable package
* parallelize build ?

### primal-dual decomposition

* implement ?
* if so, move intersecting code with dual decomposition to decomposition_tools.jl

### sparsity 

* Term Sparsity
* tests for Correlative Sparsity

### general

* normalize problem for numerical stability ?
* create Type to dispatch relaxation models and share code ?
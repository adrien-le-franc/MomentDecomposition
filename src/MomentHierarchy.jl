module MomentHierarchy

using JuMP
using DynamicPolynomials
using Combinatorics
using LinearAlgebra

# certification
using SparseArrays

# sparsity
using Graphs
using ChordalGraph

# decomposition
using Distributed
 
include("basics/utils.jl")
include("basics/polynomials.jl")
include("basics/moments.jl")

include("models/dense.jl")
include("models/decomposition_tools.jl")
include("models/correlative_sparsity.jl")
include("models/other_models/nlp.jl")
include("models/other_models/dummy_decomposition.jl")

include("sparsity/correlative.jl")

include("models/dual_decomposition/subproblems.jl")
include("models/dual_decomposition/master.jl")
include("models/dual_decomposition/oracle.jl")

include("models/certification/build_model.jl")
include("models/certification/oracle.jl")

end 

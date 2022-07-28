module MomentHierarchy

using JuMP
using DynamicPolynomials
using Combinatorics

# sparsity
using Graphs
using ChordalGraph

# decomposition
using Distributed
using SubgradientMethods # unregistered package !
const SGM = SubgradientMethods

# check utility
using LinearAlgebra
 
include("basics/utils.jl")
include("basics/polynomials.jl")
include("basics/moments.jl")

include("models/dense.jl")

include("sparsity/correlative.jl")

include("models/decomposition_tools.jl")
include("models/dual_decomposition/subproblems.jl")
include("models/dual_decomposition/master.jl")

end 

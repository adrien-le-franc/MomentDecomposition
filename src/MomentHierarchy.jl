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

include("models/hierarchies/utils.jl")
include("models/hierarchies/moment.jl")
include("models/hierarchies/sos.jl")

include("models/other_models/nlp.jl")
include("models/other_models/dummy_decomposition.jl")

include("sparsity/correlative.jl")

include("models/dual_decomposition/utils.jl")
include("models/dual_decomposition/subproblems.jl")
include("models/dual_decomposition/master.jl")
include("models/dual_decomposition/oracle.jl")

include("models/certification/build_model.jl")
include("models/certification/oracle.jl")

end 

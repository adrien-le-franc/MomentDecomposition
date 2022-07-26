module MomentDecomposition

using JuMP
using DynamicPolynomials
using Combinatorics
using LinearAlgebra
using Graphs # ??
using ChordalGraph
using Distributed
using SubgradientMethods

const SubgradientMethods = SM

include("utils.jl")
include("polynomials.jl")
include("moments.jl")
include("models/dense.jl")
include("models/decomposed.jl")
include("sparsity/correlative.jl")

end 

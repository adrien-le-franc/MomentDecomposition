# julia 1.6
#
# tools for all sparsity patterns


struct SparsityPattern{T <: Integer} 
	variable_sets::Vector{Vector{T}}
	monomial_sets::Union{Nothing, Vector{Vector{Vector{UInt16}}}}
	minimal_multipliers::Bool
end
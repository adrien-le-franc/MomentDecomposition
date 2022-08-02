# julia 1.6
#
# decomposition tools for the hierarchy of moment relaxations


function assign_constraint_to_set(polynomial::SparsePolynomial, 
	variable_sets::Vector{Vector{T}}) where T<:Integer

	for (k, set) in enumerate(variable_sets)
		if issubset(unique(vcat(polynomial.support...)), set)
			return k
		end
	end

	error("could not decompose constraint polynomial $(polynomial) over variable sets")

	return nothing

end

function pairs(variable_sets::Vector{Vector{T}}) where T<:Integer
	return enumerate(combinations(1:length(variable_sets), 2))
end

function intersect(variable_sets::Vector{Vector{T}}, 
	pair::Vector{Int64}) where T<:Integer
	return Base.intersect(variable_sets[pair[1]], variable_sets[pair[2]])
end
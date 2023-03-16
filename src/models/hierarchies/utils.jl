# julia 1.6
#
# useful functions for the moment-SOS hierarchy 


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

function update!(dict::Dict{Int64, Float64}, index::Int64, coefficient::Float64)

	if index in keys(dict)
		dict[index] += coefficient
	else
		dict[index] = coefficient
	end

	return nothing
end

function dense_relaxation_set(pop::POP)

	set = unique(vcat(pop.objective.support...))

	if pop.inequality_constraints != nothing
		for polynomial in pop.inequality_constraints
			append!(set, unique(vcat(polynomial.support...)))
		end
	end

	if pop.equality_constraints != nothing
		for polynomial in pop.equality_constraints
			append!(set, unique(vcat(polynomial.support...)))
		end
	end

	return [unique(set)]

end
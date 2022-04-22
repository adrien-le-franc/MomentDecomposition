# julia 1.6
#
# tools for moments of polynomials


# decomposed -> UInt64 ?

function moment_columns(pop::POP, max_order::Int64)
	return enumerate(with_replacement_combinations(0:pop.n_variables, max_order))
end

function moment_rows(pop::POP, max_order::Int64, column::Int64)
	return enumerate(first(with_replacement_combinations(0:pop.n_variables, max_order), column))
end

function monomial_product(alpha::Vector{Int64}, beta::Vector{Int64})
	support = vcat(alpha, beta)
	return sort!(convert.(UInt16, support[ support .!= 0]))
end

function monomial_product(alpha::Vector{Int64}, beta::Vector{Int64}, gamma::Vector{UInt16})
	support = vcat(alpha, beta)
	return sort!(vcat(convert.(UInt16, support[ support .!= 0]), gamma))
end

# decomposed -> all UInt16

function moment_columns(variables::Vector{UInt16}, max_order::Int64)
	return enumerate(with_replacement_combinations(vcat([0x0000], variables), max_order))
end

function moment_rows(variables::Vector{UInt16}, max_order::Int64, column::Int64)
	return enumerate(first(with_replacement_combinations(vcat([0x0000], variables), max_order), column))
end

function monomial_product(alpha::Vector{UInt16}, beta::Vector{UInt16})
	support = vcat(alpha, beta)
	return sort!(support[ support .!= 0x0000])
end

function monomial_product(alpha::Vector{UInt16}, beta::Vector{UInt16}, gamma::Vector{UInt16})
	support = vcat(alpha, beta)
	return sort!(vcat(support[ support .!= 0x0000], gamma))
end

function moments(variables::Vector{UInt16}, max_order::Int64)
	return with_replacement_combinations(vcat([0x0000], variables), max_order)
end

function monomial(alpha::Vector{UInt16})
	return sort!(alpha[ alpha .!= 0x0000])
end





function n_moments(variables::Vector{Int64}, relaxation_order::Int64)
	return binomial(length(variables) + relaxation_order, length(variables))
end

function n_moments(n_variables::Int64, relaxation_order::Int64)
	return binomial(n_variables + relaxation_order, n_variables)
end

function localizing_matrix_order(relaxation_order::Int64, f::SparsePolynomial)
	return relaxation_order - ceil(Int64, degree(f)/2.0)
end
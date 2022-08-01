# julia 1.6
#
# tools for moments and monomials


function n_moments(n_variables::Int64, relaxation_order::Int64)
	return binomial(n_variables + relaxation_order, n_variables)
end

function n_moments(variables::Vector{T}, relaxation_order::Int64) where T<:Integer
	return n_moments(length(variables), relaxation_order)
end

function n_moments(pop::POP, relaxation_order::Int64)
	return n_moments(pop.n_variables, relaxation_order)
end

function localizing_matrix_order(f::SparsePolynomial, relaxation_order::Int64)
	return relaxation_order - ceil(Int64, degree(f)/2.0)
end

# columns of moment/localization matrix

function moment_columns(variables::Vector{UInt16}, max_order::Int64)
	return enumerate(with_replacement_combinations(vcat([0x0000], variables), max_order))
end

function moment_columns(variables::Vector{Int64}, max_order::Int64)
	return moment_columns(convert.(UInt16, variables), max_order)
end

function moment_columns(pop::POP, max_order::Int64)
	return moment_columns(collect(0x0001:convert(UInt16, pop.n_variables)), max_order)
end

# rows of moment/localization matrix

function moment_rows(variables::Vector{UInt16}, max_order::Int64, column::Int64)
	return enumerate(first(with_replacement_combinations(vcat([0x0000], variables), max_order), column))
end

function moment_rows(variables::Vector{Int64}, max_order::Int64, column::Int64)
	return moment_rows(convert.(UInt16, variables), max_order, column)
end

function moment_rows(pop::POP, max_order::Int64, column::Int64)
	return moment_rows(collect(0x0001:convert(UInt16, pop.n_variables)), max_order, column)
end

# moments for coupling constraints

function coupling_moments(variables::Vector{UInt16}, max_order::Int64)
	return with_replacement_combinations(vcat([0x0000], variables), max_order)
end

function coupling_moments(variables::Vector{Int64}, max_order::Int64)
	return coupling_moments(convert.(UInt16, variables), max_order)
end

# monomials 

function monomial(alpha::Vector{UInt16})
	return sort!(alpha[alpha .!= 0x0000])
end

function monomial_product(alpha::Vector{UInt16}, beta::Vector{UInt16}, 
	gamma::Vector{UInt16}=UInt16[])
	return monomial(vcat(alpha, beta, gamma))
end
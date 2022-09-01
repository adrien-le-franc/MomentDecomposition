# julia 1.6
#
# dual decomposition subroblems for the hierarchy of moment relaxations


mutable struct SubProblem

	set_id::Int64
	model::Model
	polynomial_objective::AffExpr
	coupling_terms::Vector{Vector{AffExpr}}
	multiplier_ids::Vector{Int64}

end

function SubProblem(k::Int64)
	return SubProblem(k, Model(), AffExpr(0.), Vector{AffExpr}[], Int64[])
end

struct MultiplierInformation{T<:Integer}
	moments::Vector{Vector{UInt16}}
	pair::Tuple{T, T}
end

mutable struct Multiplier
	value::Vector{Vector{Float64}}
	information::Vector{MultiplierInformation}
end

function set_moment_variables!(model::Model, set::Vector{T}, 
	relaxation_order::Int64) where T<:Integer

	if !all(set .> 0)
		error("variable indices in $(set) must be strictly positive")
	end
	
	n_moment_variables = n_moments(set, 2*relaxation_order)
	@variable(model, y[1:n_moment_variables])

	return nothing

end

function set_moment_matrix!(model::Model, set::Vector{T}, 
	relaxation_order::Int64) where T<:Integer

	moment_labels = Dict{Vector{UInt16}, Int64}()
	n_monomials = n_moments(set, relaxation_order)
	M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

	count = 0

	for (j, alpha) in moment_columns(set, relaxation_order)
		for (i, beta) in moment_rows(set, relaxation_order, j)

			label = monomial_product(alpha, beta)

			if !(label in keys(moment_labels))
				count += 1
				moment_labels[label] = count
			end

			M[i, j] = model[:y][moment_labels[label]]

		end
	end

	@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

	return moment_labels

end

function set_polynomial_objective!(subproblems::Vector{SubProblem}, pop::POP, 
	moment_labels::Dict{Int64, Dict{Vector{UInt16}, Int64}}, 
	variable_sets::Vector{Vector{T}}) where T<:Integer
	
	n_sets = length(variable_sets)

	for (support, coefficient) in terms(pop.objective)
		for (k, set) in enumerate(variable_sets)

			if issubset(unique(support), set)
				add_to_expression!(subproblems[k].polynomial_objective, 
					coefficient*subproblems[k].model[:y][moment_labels[k][support]])
				break
			elseif k == n_sets
				error("could not decompose objective over variable sets")
			end

		end
	end

	return nothing

end

function set_inequality_constraint!(subproblem::SubProblem, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64}, 
	variable_sets::Vector{T}) where T<:Integer

	localizing_order = localizing_matrix_order(polynomial, relaxation_order)

	if localizing_order == 0

		@constraint(subproblem.model, 
			linear_expression(subproblem.model, polynomial, moment_labels) >= 0)

	else

		n_monomials = n_moments(variable_sets, localizing_order)
		M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

		for (j, alpha) in moment_columns(variable_sets, localizing_order)
			for (i, beta) in moment_rows(variable_sets, localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
				M[i, j] = linear_expression(subproblem.model, polynomial, labels, moment_labels)

			end
		end

		@constraint(subproblem.model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

	end

	return nothing

end

function set_equality_constraint!(subproblem::SubProblem, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_label::Dict{Vector{UInt16}, Int64}, 
	variable_sets::Vector{T}) where T<:Integer

	localizing_order = localizing_matrix_order(polynomial, relaxation_order)
	
	if localizing_order == 0

		@constraint(subproblem.model, 
			linear_expression(subproblem.model, polynomial, moment_label) == 0)

	else

		for (j, alpha) in moment_columns(variable_sets, localizing_order)
			for (i, beta) in moment_rows(variable_sets, localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) 
					for gamma in polynomial.support]
				@constraint(subproblem.model, linear_expression(subproblem.model, 
					polynomial, labels, moment_label) == 0.)

			end
		end			

	end	
	
	return nothing

end

function set_polynomial_constraints!(subproblems::Vector{SubProblem}, pop::POP, 
	relaxation_order::Int64, 
	moment_labels::Dict{Int64, Dict{Vector{UInt16}, Int64}}, 
	variable_sets::Vector{Vector{T}}) where T<:Integer

	if pop.inequality_constraints != nothing
		for polynomial in pop.inequality_constraints

			set_ids = assign_constraint_to_sets(polynomial, variable_sets)

			for k in set_ids
				set_inequality_constraint!(subproblems[k], polynomial, 
					relaxation_order, moment_labels[k], variable_sets[k])
			end

		end	
	end

	if pop.equality_constraints != nothing
		for polynomial in pop.equality_constraints

			set_ids = assign_constraint_to_sets(polynomial, variable_sets)
			
			for k in set_ids
				set_equality_constraint!(subproblems[k], polynomial, 
					relaxation_order, moment_labels[k], variable_sets[k])
			end

		end	
	end
	
	return nothing

end

function set_coupling_terms!(subproblems::Vector{SubProblem}, pair::Vector{Int64}, id::Int64,
	intersection::Vector{T}, 
	relaxation_order::Int64, 
	moment_labels::Dict{Int64, Dict{Vector{UInt16}, Int64}}, 
	max_coupling_order::Int64) where T<:Integer

	k_1, k_2 = pair
	pair_information = MultiplierInformation(Vector{UInt16}[], (k_1, k_2))

	coupling_terms_k_1 = AffExpr[] 
	coupling_terms_k_2 = AffExpr[]

	for alpha in coupling_moments(intersection, 2*relaxation_order)

		label = monomial(alpha)

		if length(label) > max_coupling_order || length(label) == 0
			continue
		end

		push!(coupling_terms_k_1, 1*subproblems[k_1].model[:y][moment_labels[k_1][label]])
		push!(coupling_terms_k_2, -1*subproblems[k_2].model[:y][moment_labels[k_2][label]])
		push!(pair_information.moments, label)

	end

	push!(subproblems[k_1].coupling_terms, coupling_terms_k_1)
	push!(subproblems[k_2].coupling_terms, coupling_terms_k_2)
	push!(subproblems[k_1].multiplier_ids, id)
	push!(subproblems[k_2].multiplier_ids, id)

	return pair_information

end

function set_Lagrange_multipliers!(subproblems::Vector{SubProblem}, 
	relaxation_order::Int64, moment_labels::Dict{Int64, Dict{Vector{UInt16}, Int64}}, 
	variable_sets::Vector{Vector{T}}, max_coupling_order::Int64) where T<:Integer

	multiplier_information = MultiplierInformation[]
	pair_id = 0

	for pair in pairs(variable_sets)

		intersection = intersect(variable_sets, pair)

		if isempty(intersection)
			continue
		end

		pair_id += 1

		pair_information = set_coupling_terms!(subproblems, pair, pair_id, intersection, 
			relaxation_order, moment_labels, max_coupling_order)

		push!(multiplier_information, pair_information)

	end

	return multiplier_information

end
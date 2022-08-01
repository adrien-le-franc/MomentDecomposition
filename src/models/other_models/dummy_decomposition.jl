# julia 1.6
#
# dummy (a single JuMP model) decomposition for the hierarchy of moment relaxation


y_(k::Int64) = Symbol("y_$k")

function set_moment_variables!(model::Model, set::Vector{T},
	k::Int64, relaxation_order::Int64) where T<:Integer

	if !all(set .> 0)
		error("variable indices in $(set) must be strictly positive")
	end
	
	n_moment_variables = n_moments(set, 2*relaxation_order)
	y = @variable(model, [1:n_moment_variables], base_name = "y_$k")
	model[y_(k)] = y

	return nothing

end

function set_moment_matrix!(model::Model, set::Vector{T}, k::Int64, 
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

			M[i, j] = model[y_(k)][moment_labels[label]]

		end
	end

	@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

	return moment_labels

end

function set_probability_measure_constraint!(model::Model, 
	moment_labels::Dict{Vector{UInt16}, Int64}, k::Int64)

	@constraint(model, model[y_(k)][moment_labels[k][UInt16[]]] == 1.)
	return nothing

end

function set_objective!(model::Model, pop::POP, moment_labels::Dict{Vector{UInt16}, Int64}, 
	variable_sets::Vector{Vector{T}}) where T<:Integer

	f = AffExpr(0.)

	for (support, coefficient) in terms(pop.objective)
		for (k, set) in enumerate(variable_sets)

			if issubset(unique(support), set)
				add_to_expression!(f, coefficient*model[y_(k)][moment_labels[k][support]])
				break
			elseif k == n_sets
				error("could not decompose objective over variable sets")
			end

		end

	end

	@objective(model, Min, f)

	return nothing

end

function linear_expression(model::Model, f::SparsePolynomial, 
    moment_labels_k::Dict{Vector{UInt16}, Int64}, k::Int64)
    
    return sum(coefficient*model[y_(k)][moment_labels_k[support]]
                for (support, coefficient) in terms(f))

end

function linear_expression(model::Model, f::SparsePolynomial, 
    labels::Vector{Vector{UInt16}}, moment_labels_k::Dict{Vector{UInt16}, Int64}, k::Int64)
    
    return sum(coefficient*model[y_(k)][moment_labels_k[labels[i]]] 
                for (i, coefficient) in enumerate(f.coefficients))

end

function set_inequality_constraints!(model::Model, pop, relaxation_order, moment_labels, variable_sets::Vector{Vector{UInt16}})

	if pop.inequality_constraints == nothing
		return nothing
	end

	for polynomial in pop.inequality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, linear_expression(model, polynomial, moment_labels[k], k) >= 0)

		else

			n_monomials = n_moments(length(variable_sets[k]), relaxation_order)
			M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					M[i, j] = linear_expression(model, polynomial, labels, moment_labels[k], k)

				end
			end

			@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

		end

	end

	return nothing

end

function set_equality_constraints!(model::Model, pop, relaxation_order, moment_labels, variable_sets)

	if pop.equality_constraints == nothing
		return nothing
	end

	for polynomial in pop.equality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, linear_expression(model, polynomial, moment_labels[k], k) == 0)

		else

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(model, linear_expression(model, polynomial, 
						labels, moment_labels[k], k) == 0.)

				end
			end			

		end


	end

	return nothing

end

function set_coupling_constraints!(model::Model, relaxation_order::Int64, 
	moment_labels, 
	variable_sets::Vector{Vector{Integer}}, 
	max_coupling_order::Int64)

	for pair in combinations(1:length(variable_sets), 2)

		k_1, k_2 = pair
		intersection = intersect(variable_sets[k_1], variable_sets[k_2])

		for alpha in coupling_moments(intersection, 2*relaxation_order)

			label = monomial(alpha)

			if length(label) > max_coupling_order || length(label) == 0
				continue
			end

			@constraint(model, model[y_(k_1)][moment_labels[k_1][label]] == model[y_(k_2)][moment_labels[k_2][label]])

		end

	end

	return nothing 

end

function dummy_decomposed_relaxation_model(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
	max_coupling_order::Int64=relaxation_order)

	model = Model() 

	n_sets = length(variable_sets)
	moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:n_sets)
	
	for (k, set) in enumerate(variable_sets)

		set_moment_variables!(model, relaxation_order, set, k)
		moment_labels[k] = set_moment_matrix!(model, relaxation_order, set, k)
		set_probability_measure_constraint!(model, moment_labels, k)

	end

	set_objective!(model, pop, moment_labels, variable_sets)
	set_inequality_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)
	set_equality_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)

	set_coupling_constraints!(model, relaxation_order, moment_labels, variable_sets, max_coupling_order)

	return model

end
# julia 1.6
#
# decomposition tools for the hierarchy of moment relaxations


## tools for defining separable subproplems


function set_moment_variables!(model::Model, relaxation_order, set::Vector{Integer}, k::Int64)

	if !all(set .> 0)
		error("variable indices in $(set) must be strictly positive")
	end
	
	n_moment_variables = n_moments(set, 2*relaxation_order)
	y = @variable(model, [1:n_moment_variables], base_name = "y_$k")
	model[y_(k)] = y

	return nothing

end

function set_moment_matrix!(model::Model, relaxation_order::Int64, set::Vector{Integer}, k::Int64)

	moment_labels = Dict{Vector{UInt16}, Int64}()

	n_monomials = n_moments(set, relaxation_order)
	moment_matrix = @variable(model, [1:n_monomials, 1:n_monomials], PSD, base_name = "moment_matrix_$k")
	model[Symbol("moment_matrix_$k")] = moment_matrix	

	count = 0

	for (j, alpha) in moment_columns(set, relaxation_order)
		for (i, beta) in moment_rows(set, relaxation_order, j)

			label = monomial_product(alpha, beta)

			if !(label in keys(moment_labels))
				count += 1
				moment_labels[label] = count
			end

			@constraint(model, moment_matrix[i, j] == model[y_(k)][moment_labels[label]])

		end
	end

	return moment_labels

end

function set_probability_measure_constraint!(model::Model, moment_labels, k)
	@constraint(model, model[y_(k)][moment_labels[k][UInt16[]]] == 1.)
	return nothing
end

function objective_decomposition(pop::POP, variable_sets)

	n_monomials = length(pop.objective.coefficients)
	n_sets = length(variable_sets)

	coefficients_per_set = [Float64[] for k in 1:n_sets]
	support_per_set = [Vector{UInt16}[] for k in 1:n_sets]

	for n in 1:n_monomials

		monomial_support = pop.objective.support[n]

		for (k, set) in enumerate(variable_sets)
			if issubset(unique(monomial_support), set)
				push!(coefficients_per_set[k], pop.objective.coefficients[n])
				push!(support_per_set[k], monomial_support)
				break
			elseif k == n_sets
				error("could not decompose objective over variable sets")
			end
		end

	end

	return [SparsePolynomial(support_per_set[k], coefficients_per_set[k]) for k in 1:n_sets]

end

function assign_constraint_to_set(polynomial, variable_sets)

	for (k, set) in enumerate(variable_sets)
		if issubset(unique(vcat(polynomial.support...)), set)
			return k
		end
	end

	error("could not decompose constraint polynomial $(polynomial) over variable sets")

	return nothing

end


## almost decomposed relaxation (encoded in a single JuMP model, mostly for test purpose)


function set_objective!(model::Model, pop, moment_labels, variable_sets)

	objective = objective_decomposition(pop, variable_sets)

	f = AffExpr(0.) 

	for k in 1:length(variable_sets)
		n_temrs = length(objective[k].coefficients)
		if n_temrs > 0
			moments_f_k = [moment_labels[k][objective[k].support[n]] for n in 1:n_temrs]
			f += sum(objective[k].coefficients[n]*model[Symbol("y_$k")][moments_f_k[n]] for n in 1:n_temrs)	
		end
	end

	@objective(model, Min, f)

	return nothing

end

function set_inequality_constraints!(model::Model, pop, relaxation_order, moment_labels, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.inequality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, sum(polynomial.coefficients[n]*model[Symbol("y_$k")][moment_labels[k][polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) >= 0)

		else

			n_monomials = n_moments(length(variable_sets[k]), relaxation_order)
			localizing_matrix = @variable(model, [1:n_monomials, 1:n_monomials], PSD)

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(model, localizing_matrix[i, j] == sum(polynomial.coefficients[n]*model[Symbol("y_$k")][moment_labels[k][labels[n]]] for n in 1:length(polynomial.coefficients)))

				end
			end

		end


	end

	return nothing

end

function set_equality_constraints!(model::Model, pop, relaxation_order, moment_labels, variable_sets)

	for polynomial in pop.equality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, sum(polynomial.coefficients[n]*model[Symbol("y_$k")][moment_labels[k][polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) == 0)

		else

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(model, sum(polynomial.coefficients[n]*model[Symbol("y_$k")][moment_labels[k][labels[n]]] for n in 1:length(polynomial.coefficients)) == 0.)

				end
			end			

		end


	end

	return nothing

end

function set_coupling_constraints!(model::Model, relaxation_order, moment_labels, 
	variable_sets::Vector{Vector{Integer}}, max_coupling_order)

	for pair in combinations(1:length(variable_sets), 2)

		k_1, k_2 = pair
		intersection = intersect(variable_sets[k_1], variable_sets[k_2])

		for alpha in coupling_moments(intersection, 2*relaxation_order) # remove moment 0...0 ?

			label = monomial(alpha)

			if length(label) > max_coupling_order
				continue
			end

			@constraint(model, model[Symbol("y_$(k_1)")][moment_labels[k_1][label]] == model[Symbol("y_$(k_2)")][moment_labels[k_2][label]])

		end

	end

	return nothing 

end

function almost_decomposed_relaxation(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
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
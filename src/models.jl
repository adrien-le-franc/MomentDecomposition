# julia 1.6
#
# models for the hierarchy of moment relaxations


function set_moment_variables!(model, pop, relaxation_order)
	
	n_moment_variables = n_moments(pop.n_variables, 2*relaxation_order)
	@variable(model, y[1:n_moment_variables])

	return nothing

end

function set_moment_matrix!(model, pop, relaxation_order)

	moment_labels = Dict{Vector{Int64}, Int64}()

	n_monomials = n_moments(pop.n_variables, relaxation_order)
	@variable(model, moment_matrix[1:n_monomials, 1:n_monomials], PSD)

	count = 0

	for (j, alpha) in moment_columns(pop, relaxation_order)
		for (i, beta) in moment_rows(pop, relaxation_order, j)

			label = monomial_product(alpha, beta)

			if !(label in keys(moment_labels))
				count += 1
				moment_labels[label] = count
			end

			@constraint(model, moment_matrix[i, j] == model[:y][moment_labels[label]])

		end
	end

	return moment_labels

end

function set_objective!(model, pop, moment_labels)

	n_temrs = length(pop.objective.support)
	moments_f = [moment_labels[pop.objective.support[n]] for n in 1:n_temrs]
	f = sum(pop.objective.coefficients[n]*model[:y][moments_f[n]] for n in 1:n_temrs)

	@objective(model, Min, f)

	return nothing

end

function set_probability_measure_constraint!(model, moment_labels)
	@constraint(model, model[:y][moment_labels[UInt16[]]] == 1.)
	return nothing
end

function set_inequality_constraints!(model, pop, relaxation_order, moment_labels)

	for polynomial in pop.inequality_constraints

		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, sum(polynomial.coefficients[k]*model[:y][moment_labels[polynomial.support[k]]]
				for k in 1:length(polynomial.coefficients)) >= 0)

		else

			n_monomials = n_moments(pop.n_variables, localizing_order)
			localizing_matrix = @variable(model, [1:n_monomials, 1:n_monomials], PSD)

			for (j, alpha) in moment_columns(pop, localizing_order)
				for (i, beta) in moment_rows(pop, localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(model, localizing_matrix[i, j] == sum(polynomial.coefficients[k]*model[:y][moment_labels[labels[k]]] for k in 1:length(polynomial.coefficients)))

				end
			end

		end


	end

	return nothing

end

function set_equality_constraints!(model, pop, relaxation_order, moment_labels)

	for polynomial in pop.equality_constraints

		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, sum(polynomial.coefficients[k]*model[:y][moment_labels[polynomial.support[k]]]
				for k in 1:length(polynomial.coefficients)) == 0)

		else

			for (j, alpha) in moment_columns(pop, localizing_order)
				for (i, beta) in moment_rows(pop, localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(model, sum(polynomial.coefficients[k]*model[:y][moment_labels[labels[k]]] for k in 1:length(polynomial.coefficients)) == 0.)

				end
			end			

		end

	end

	return nothing

end

function dense_relaxation(pop::POP, relaxation_order::Int64)

	model = Model() # struct for moment hierarchy ?

	set_moment_variables!(model, pop, relaxation_order)
	moment_labels = set_moment_matrix!(model, pop, relaxation_order)
	set_objective!(model, pop, moment_labels)
	set_probability_measure_constraint!(model, moment_labels)

	set_inequality_constraints!(model, pop, relaxation_order, moment_labels)
	set_equality_constraints!(model, pop, relaxation_order, moment_labels)

	return model

end






















function set_moment_variables!(model, relaxation_order, set::Vector{UInt16}, k::Int64)
	
	n_moment_variables = n_moments(length(set), 2*relaxation_order)
	y = @variable(model, [1:n_moment_variables], base_name = "y_$k")
	model[Symbol("y_$k")] = y

	return nothing

end

function set_moment_matrix!(model::JuMP.Model, relaxation_order::Int64, set::Vector{UInt16}, k::Int64)

	moment_labels = Dict{Vector{UInt16}, Int64}()

	n_monomials = n_moments(length(set), relaxation_order)
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

			@constraint(model, moment_matrix[i, j] == model[Symbol("y_$k")][moment_labels[label]])

		end
	end

	return moment_labels

end

function set_probability_measure_constraint!(model, moment_labels, set, k)
	@constraint(model, model[Symbol("y_$k")][moment_labels[k][UInt16[]]] == 1.)
	return nothing
end

function set_objective!(model, pop, moment_labels, variable_sets)

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

	f = AffExpr(0.) 

	for k in 1:n_sets
		n_temrs = length(coefficients_per_set[k])
		if n_temrs > 0
			moments_f_k = [moment_labels[k][support_per_set[k][n]] for n in 1:n_temrs]
			f += sum(coefficients_per_set[k][n]*model[Symbol("y_$k")][moments_f_k[n]] for n in 1:n_temrs)	
		end
	end

	@objective(model, Min, f)

	return nothing

end

function set_inequality_constraints!(model, pop, relaxation_order, moment_labels, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.inequality_constraints

		# assign constraint to set
		k = 0
		for (n, set) in enumerate(variable_sets)
			if issubset(unique(vcat(polynomial.support...)), set)
				k = n
				break
			end
		end
		if k == 0
			error("could not decompose $(polynomial) >= 0 over variable sets")
		end

		# write constraint
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

function set_equality_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)

	for polynomial in pop.equality_constraints

		# assign constraint to set
		k = 0
		for (n, set) in enumerate(variable_sets)
			if issubset(unique(vcat(polynomial.support...)), set)
				k = n
				break
			end
		end
		if k == 0
			error("could not decompose $(polynomial) >= 0 over variable sets")
		end

		# write constraint
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

function set_coupling_constraints!(model, relaxation_order, moment_labels, variable_sets)

	for pair in combinations(1:length(variable_sets), 2)

		k_1, k_2 = pair
		intersection = intersect(variable_sets[k_1], variable_sets[k_2])

		for alpha in moments(intersection, 2*relaxation_order)

			label = monomial(alpha)
			@constraint(model, model[Symbol("y_$(k_1)")][moment_labels[k_1][label]] == model[Symbol("y_$(k_2)")][moment_labels[k_2][label]])

		end

	end

	return nothing 

end

function decomposed_relaxation(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{Int64}})

	model = Model() 

	n_sets = length(variable_sets)
	moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:n_sets)
	uint16_variable_sets = [convert.(UInt16, set) for set in variable_sets] 

	for (k, set) in enumerate(uint16_variable_sets)

		set_moment_variables!(model, relaxation_order, set, k)
		moment_labels[k] = set_moment_matrix!(model, relaxation_order, set, k)
		set_probability_measure_constraint!(model, moment_labels, set, k)

	end

	set_objective!(model, pop, moment_labels, uint16_variable_sets)
	set_inequality_constraints!(model, pop, relaxation_order, moment_labels, uint16_variable_sets)
	set_equality_constraints!(model, pop, relaxation_order, moment_labels, uint16_variable_sets)

	set_coupling_constraints!(model, relaxation_order, moment_labels, uint16_variable_sets)

	return model

end
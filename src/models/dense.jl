# julia 1.6
#
# dense models for the hierarchy of moment relaxations


# dense relaxation

function set_moment_variables!(model, pop, relaxation_order)
	
	n_moment_variables = n_moments(pop, 2*relaxation_order)
	@variable(model, y[1:n_moment_variables])

	return nothing

end

function set_moment_matrix!(model, pop, relaxation_order)

	moment_labels = Dict{Vector{UInt16}, Int64}() 

	n_monomials = n_moments(pop, relaxation_order)

	M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

	count = 0

	for (j, alpha) in moment_columns(pop, relaxation_order)
		for (i, beta) in moment_rows(pop, relaxation_order, j)

			label = monomial_product(alpha, beta)

			if !(label in keys(moment_labels))
				count += 1
				moment_labels[label] = count
			end

			M[i, j] = AffExpr(0.)
			add_to_expression!(M[i, j], model[:y][moment_labels[label]])

		end
	end

	@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

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

	if pop.inequality_constraints == nothing
		return nothing
	end

	for polynomial in pop.inequality_constraints

		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(model, sum(polynomial.coefficients[k]*model[:y][moment_labels[polynomial.support[k]]]
				for k in 1:length(polynomial.coefficients)) >= 0)

		else

			n_monomials = n_moments(pop.n_variables, localizing_order)

			M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

			for (j, alpha) in moment_columns(pop, localizing_order)
				for (i, beta) in moment_rows(pop, localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					

					M[i, j] = AffExpr(0.)
					add_to_expression!(M[i, j], sum(polynomial.coefficients[k]*model[:y][moment_labels[labels[k]]] for k in 1:length(polynomial.coefficients)))

				end
			end


			@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())


		end


	end

	return nothing

end

function set_equality_constraints!(model, pop, relaxation_order, moment_labels)

	if pop.equality_constraints == nothing
		return nothing
	end

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













# NLP


function polynomial_expression(model, polynomial::SparsePolynomial)

	support = [convert.(Int64, monomial) for monomial in polynomial.support]
	n_terms = length(support)

	expression = @NLexpression(model, 
		sum(polynomial.coefficients[i]*prod(model[:x][k] for k in support[i]) for i in 1:n_terms))

	return expression

end

function non_linear_problem(pop::POP)

	model = Model()

	# var
	@variable(model, x[1:pop.n_variables])

	# cons

	for g_inequality in pop.inequality_constraints
		polynomial = polynomial_expression(model, g_inequality)
		@NLconstraint(model, polynomial >= 0.)
	end

	for g_equality in pop.equality_constraints
		polynomial = polynomial_expression(model, g_equality)
		@NLconstraint(model, polynomial == 0.)
	end

	# obj
	polynomial = polynomial_expression(model, pop.objective)
	@NLobjective(model, Min, polynomial)
	
	return model

end
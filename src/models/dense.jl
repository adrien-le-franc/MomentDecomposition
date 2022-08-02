# julia 1.6
#
# dense models for the hierarchy of moment relaxations


function set_moment_variables!(model::Model, pop::POP, relaxation_order::Int64)
	
	n_moment_variables = n_moments(pop, 2*relaxation_order)
	@variable(model, y[1:n_moment_variables])

	return nothing

end

function set_moment_matrix!(model::Model, pop::POP, relaxation_order::Int64)

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

			M[i, j] = model[:y][moment_labels[label]]

		end
	end

	@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

	return moment_labels

end

function linear_expression(model::Model, f::SparsePolynomial, 
    moment_labels::Dict{Vector{UInt16}, Int64})
    
    return sum(coefficient*model[:y][moment_labels[support]]
                for (support, coefficient) in terms(f))

end

function linear_expression(model::Model, f::SparsePolynomial, 
    labels::Vector{Vector{UInt16}}, moment_labels::Dict{Vector{UInt16}, Int64})
    
    return sum(coefficient*model[:y][moment_labels[labels[i]]] 
                for (i, coefficient) in enumerate(f.coefficients))

end

function set_objective!(model::Model, pop::POP, moment_labels::Dict{Vector{UInt16}, Int64})

	@objective(model, Min, linear_expression(model, pop.objective, moment_labels))
	return nothing

end

function set_probability_measure_constraint!(model::Model, 
	moment_labels::Dict{Vector{UInt16}, Int64})
	
	@constraint(model, model[:y][moment_labels[UInt16[]]] == 1.)
	return nothing

end

function set_inequality_constraint!(model::Model, pop::POP, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64})
	
	localizing_order = localizing_matrix_order(polynomial, relaxation_order)

	if localizing_order == 0

		@constraint(model, linear_expression(model, polynomial, moment_labels) >= 0)

	else

		n_monomials = n_moments(pop.n_variables, localizing_order)
		M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

		for (j, alpha) in moment_columns(pop, localizing_order)
			for (i, beta) in moment_rows(pop, localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
				M[i, j] = linear_expression(model, polynomial, labels, moment_labels)

			end
		end

		@constraint(model, LinearAlgebra.Symmetric(M) >= 0, PSDCone())

	end

	return nothing

end

function set_equality_constraint!(model::Model, pop::POP, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64})
	
	localizing_order = localizing_matrix_order(polynomial, relaxation_order)

	if localizing_order == 0

		@constraint(model, linear_expression(model, polynomial, moment_labels) == 0)

	else

		for (j, alpha) in moment_columns(pop, localizing_order)
			for (i, beta) in moment_rows(pop, localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
				@constraint(model, linear_expression(model, polynomial, 
					labels, moment_labels) == 0.)

			end
		end			

	end

	return nothing

end

function set_polynomial_constraints!(model::Model, pop::POP, relaxation_order::Int64, 
	moment_labels::Dict{Vector{UInt16}, Int64})

	if pop.inequality_constraints != nothing
		for polynomial in pop.inequality_constraints
			set_inequality_constraint!(model, pop, polynomial, relaxation_order, moment_labels)

		end	
	end

	if pop.equality_constraints != nothing
		for polynomial in pop.equality_constraints
			set_equality_constraint!(model, pop, polynomial, relaxation_order, moment_labels)

		end	
	end
	

	return nothing

end

function dense_relaxation_model(pop::POP, relaxation_order::Int64)

	model = Model()

	set_moment_variables!(model, pop, relaxation_order)
	moment_labels = set_moment_matrix!(model, pop, relaxation_order)
	set_objective!(model, pop, moment_labels)
	set_probability_measure_constraint!(model, moment_labels)
	set_polynomial_constraints!(model, pop, relaxation_order, moment_labels)

	return model

end
# julia 1.6
#
# generic sparse model for the hierarchy of moment relaxations


function set_moment_matrices!(model::Model, pop::POP, relaxation_order::Int64, 
	sets::Vector{Vector{T}}) where T <: Integer

	moment_labels = Dict{Vector{UInt16}, Int64}() 
	n_monomials = n_moments(pop, relaxation_order)
	y = Vector{VariableRef}()
	count = 0

	for set in sets

		n_monomials = n_moments(set, relaxation_order)
		M = Matrix{AffExpr}(undef, n_monomials, n_monomials)
		
		for (j, alpha) in moment_columns(set, relaxation_order)
			for (i, beta) in moment_rows(set, relaxation_order, j)

				label = monomial_product(alpha, beta)

				if !(label in keys(moment_labels))

					count += 1
					push!(y, @variable(model)) #, base_name="y_$(count)")) # base_name can be costly for large models: remove ?
					moment_labels[label] = count
				
				end

				M[i, j] = y[moment_labels[label]]

			end
		end

		@constraint(model, Symmetric(M) >= 0, PSDCone())

	end

	model[:y] = y

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

function set_inequality_constraint!(model::Model, pop::POP, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64}, 
	variable_sets::Vector{Vector{T}}) where T <: Integer

	localizing_order = localizing_matrix_order(polynomial, relaxation_order)

	if localizing_order == 0

		@constraint(model, linear_expression(model, polynomial, moment_labels) >= 0)

	else

		k = assign_constraint_to_set(polynomial, variable_sets)
		n_monomials = n_moments(variable_sets[k], localizing_order)
		M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

		for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
			for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
				M[i, j] = linear_expression(model, polynomial, labels, moment_labels)

			end
		end

		@constraint(model, Symmetric(M) >= 0, PSDCone())

	end

	return nothing

end

function set_equality_constraint!(model::Model, pop::POP, polynomial::SparsePolynomial, 
	relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64},
	variable_sets::Vector{Vector{T}}) where T <: Integer
	
	localizing_order = localizing_matrix_order(polynomial, relaxation_order)

	if localizing_order == 0

		@constraint(model, linear_expression(model, polynomial, moment_labels) == 0)

	else

		k = assign_constraint_to_set(polynomial, variable_sets)

		for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
			for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

				labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
				@constraint(model, linear_expression(model, polynomial, 
					labels, moment_labels) == 0.)

			end
		end			

	end

	return nothing

end

function set_polynomial_constraints!(model::Model, pop::POP, relaxation_order::Int64, 
	moment_labels::Dict{Vector{UInt16}, Int64},
	variable_sets::Vector{Vector{T}}) where T <: Integer

	if pop.inequality_constraints != nothing
		for polynomial in pop.inequality_constraints
			set_inequality_constraint!(model, pop, polynomial, relaxation_order, 
				moment_labels, variable_sets)
		end	
	end

	if pop.equality_constraints != nothing
		for polynomial in pop.equality_constraints
			set_equality_constraint!(model, pop, polynomial, relaxation_order, 
				moment_labels, variable_sets)
		end	
	end
	
	return nothing

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

function moment_relaxation_model(pop::POP, relaxation_order::Int64, 
	sparsity_sets::Vector{Vector{T}}=dense_relaxation_set(pop)) where T <: Integer

	model = Model()

	moment_labels = set_moment_matrices!(model, pop, relaxation_order, sparsity_sets)
	set_objective!(model, pop, moment_labels)
	set_probability_measure_constraint!(model, moment_labels)
	set_polynomial_constraints!(model, pop, relaxation_order, moment_labels, sparsity_sets)

	return model

end
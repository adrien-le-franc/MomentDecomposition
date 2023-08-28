# julia 1.6
#
# generic sparse model for the hierarchy of moment relaxations


function set_moment_matrices!(model::Model, variable_sets::Vector{Vector{T}}, monomial_sets, relaxation_order::Int64) where T <: Integer

	moment_labels = Dict{Vector{UInt16}, Int64}() 
	y = Vector{VariableRef}()
	count = 0

	for (k, set) in enumerate(variable_sets)

		blocks = compute_tsp_blocks_X_0(set, monomial_sets[k], relaxation_order)

		for (l, block) in matrix_blocks(blocks)

			M = Matrix{AffExpr}(undef, length(block), length(block))	

			for (j, (_, alpha)) in moment_columns(block, set, relaxation_order)
				for (i, (_, beta)) in moment_rows(block, set, relaxation_order, j)

					monomial = monomial_product(alpha, beta)

					if !(monomial in keys(moment_labels))
						count += 1
						moment_labels[monomial] = count
						push!(y, @variable(model)) 
					end

					M[i, j] = y[moment_labels[monomial]]					

				end
			end

			@constraint(model, Symmetric(M) >= 0, PSDCone())

		end

		for (l, (_, alpha)) in moment_columns(diagonal_block(blocks), set, relaxation_order)

			monomial = monomial_product(alpha, alpha)

			if !(monomial in keys(moment_labels))
				count += 1
				moment_labels[monomial] = count
				push!(y, @variable(model)) 
			end

			@constraint(model, y[moment_labels[monomial]] >= 0.)

		end

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
	count = length(moment_labels)

	if localizing_order == 0

		expression = AffExpr(0.)

		for (support, coefficient) in terms(polynomial)

			if !(support in keys(moment_labels))
				count += 1
				moment_labels[support] = count
				push!(model[:y], @variable(model))
			end

			add_to_expression!(expression, coefficient, model[:y][moment_labels[support]])

		end

		@constraint(model, expression >= 0)

	else

		#k = assign_constraint_to_set(polynomial, variable_sets)
		#set = unique(vcat(polynomial.support...))
		set = unique(vcat(polynomial.support...))

		n_monomials = n_moments(set, localizing_order)
		M = Matrix{AffExpr}(undef, n_monomials, n_monomials)

		for (j, alpha) in moment_columns(set, localizing_order)
			for (i, beta) in moment_rows(set, localizing_order, j)

				expression = AffExpr(0.)

				for (gamma, coefficient) in terms(polynomial)

					monomial = monomial_product(alpha, beta, gamma)

					if !(monomial in keys(moment_labels))
						count += 1
						moment_labels[monomial] = count
						push!(model[:y], @variable(model))
					end

					add_to_expression!(expression, coefficient, model[:y][moment_labels[monomial]])

				end

				M[i, j] = expression

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
	count = length(moment_labels)

	if localizing_order == 0

		expression = AffExpr(0.)

		for (support, coefficient) in terms(polynomial)

			if !(support in keys(moment_labels))
				count += 1
				moment_labels[support] = count
				push!(model[:y], @variable(model))
			end

			add_to_expression!(expression, coefficient, model[:y][moment_labels[support]])

		end

		@constraint(model, expression == 0)

	else

		#k = assign_constraint_to_set(polynomial, variable_sets)
		#set = unique(vcat(polynomial.support...))
		set = unique(vcat(polynomial.support...))

		for (j, alpha) in moment_columns(set, localizing_order)
			for (i, beta) in moment_rows(set, localizing_order, j)

				expression = AffExpr(0.)

				for (gamma, coefficient) in terms(polynomial)

					monomial = monomial_product(alpha, beta, gamma)

					if !(monomial in keys(moment_labels))
						count += 1
						moment_labels[monomial] = count
						push!(model[:y], @variable(model))
					end

					add_to_expression!(expression, coefficient, model[:y][moment_labels[monomial]])

				end

				@constraint(model, expression == 0)

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
	variable_sets::Vector{Vector{T}},  
	# =dense_relaxation_set(pop),
	monomial_sets::Vector{Vector{Vector{UInt16}}};
	minimal_multipliers=true,
	moment_scaling=false,
	return_monomials=false) where T <: Integer

	model = Model()

	moment_labels = set_moment_matrices!(model, variable_sets, monomial_sets, relaxation_order)
	set_objective!(model, pop, moment_labels)
	set_probability_measure_constraint!(model, moment_labels)
	set_polynomial_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)

	return model

end
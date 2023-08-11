# julia 1.6.5


function set_X_0!(model, variable_sets, monomial_sets, relaxation_order)

	monomial_index = Dict{Vector{UInt16}, Int64}()
	X_0 = Vector{Vector{Symmetric{VariableRef, Matrix{VariableRef}}}}(undef, length(variable_sets))
	x_0 = Vector{Vector{VariableRef}}(undef, length(variable_sets))
	linear_operator = Vector{AffExpr}()
	
	count = 0

	for (k, set) in enumerate(variable_sets)

		# changes !!

		blocks = compute_tsp_blocks_X_0(set, monomial_sets[k], relaxation_order)

		n_matrices = sum([length.(blocks) .> 1]...)
		X_0_k = Vector{Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n_matrices)
		
		for (l, block) in matrix_blocks(blocks)

			X_0_k[l] = @variable(model, [1:length(block), 1:length(block)], PSD)

			for (j, (_, alpha)) in moment_columns(block, set, relaxation_order)
				for (i, (_, beta)) in moment_rows(block, set, relaxation_order, j)

					monomial = monomial_product(alpha, beta)

					if !(monomial in keys(monomial_index))
						count += 1
						monomial_index[monomial] = count
						push!(linear_operator, AffExpr(0., X_0_k[l][i, j] => 2 - (i==j))) 
					else
						add_to_expression!(linear_operator[monomial_index[monomial]], 2 - (i==j), X_0_k[l][i, j])
					end

				end
			end

		end

		X_0[k] = X_0_k

		n_scalars = length(blocks) - n_matrices		
		x_0_k =	Vector{VariableRef}(undef, n_scalars)

		for (l, (_, alpha)) in moment_columns(diagonal_block(blocks), set, relaxation_order)

			x_0_k[l] =  @variable(model)
			@constraint(model, x_0_k[l] >= 0.)

			monomial = monomial_product(alpha, alpha)

			if !(monomial in keys(monomial_index))
				count += 1
				monomial_index[monomial] = count
				push!(linear_operator, AffExpr(0., x_0_k[l] => 1)) 
			else
				add_to_expression!(linear_operator[monomial_index[monomial]], 1, x_0_k[l])
			end

		end

		x_0[k] = x_0_k

	end

	model[:X_0] = X_0
	model[:x_0] = x_0
	model[:linear_operator] = linear_operator

	return monomial_index

end

function set_X_j!(model, pop, monomial_index, sets, relaxation_order)

	count = length(monomial_index) 

	# inequality constraints

	if pop.inequality_constraints != nothing

		X_ineq = Vector{Union{VariableRef, Symmetric{VariableRef, Matrix{VariableRef}}}}(undef, length(pop.inequality_constraints))

		for (l, polynomial) in enumerate(pop.inequality_constraints)

			localizing_order = localizing_matrix_order(polynomial, relaxation_order)
			k = assign_constraint_to_set(polynomial, sets)
			set = sets[k]

			if localizing_order == 0

				X_ineq[l] =  @variable(model)
				@constraint(model, X_ineq[l] >= 0.)

				for (support, coefficient) in terms(polynomial)

					if !(support in keys(monomial_index))
						count += 1
						monomial_index[support] = count
						push!(model[:linear_operator], AffExpr(0., X_ineq[l] => coefficient)) 
					else
						add_to_expression!(model[:linear_operator][monomial_index[support]], coefficient, X_ineq[l])
					end

				end

			else

				matrix_size = n_moments(set, localizing_order)
				X_ineq[l] = @variable(model, [1:matrix_size, 1:matrix_size], PSD)

				for (j, alpha) in moment_columns(set, localizing_order)
					for (i, beta) in moment_rows(set, localizing_order, j)

						monomials = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support] 

						for (n, monomial) in enumerate(monomials)

							if !(monomial in keys(monomial_index))
								count += 1
								monomial_index[monomial] = count
								push!(model[:linear_operator], AffExpr(0., X_ineq[l][i, j] => (2 - (i==j))*polynomial.coefficients[n])) 
							else
								add_to_expression!(model[:linear_operator][monomial_index[monomial]], (2 - (i==j))*polynomial.coefficients[n], X_ineq[l][i, j])
							end							

						end
						
					end
				end

			end

		end

		model[:X_ineq] = X_ineq
		
	end

	# equality constraints

	if pop.equality_constraints != nothing

		X_eq = Vector{Union{VariableRef, Symmetric{VariableRef, Matrix{VariableRef}}}}(undef, length(pop.equality_constraints))

		for (l, polynomial) in enumerate(pop.equality_constraints)

			localizing_order = localizing_matrix_order(polynomial, relaxation_order)
			k = assign_constraint_to_set(polynomial, sets)
			set = sets[k]

			if localizing_order == 0

				X_eq[l] = @variable(model)

				for (support, coefficient) in terms(polynomial)
					
					if !(support in keys(monomial_index))
						count += 1
						monomial_index[support] = count
						push!(model[:linear_operator], AffExpr(0., X_eq[l] => coefficient)) 
					else
						add_to_expression!(model[:linear_operator][monomial_index[support]], coefficient, X_eq[l])
					end
	
				end

			else

				matrix_size = n_moments(set, localizing_order)
				X_eq[l] = @variable(model, [1:matrix_size, 1:matrix_size], Symmetric)

				for (j, alpha) in moment_columns(set, localizing_order)
					for (i, beta) in moment_rows(set, localizing_order, j)

						monomials = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support] 

						for (n, monomial) in enumerate(monomials)
							
							if !(monomial in keys(monomial_index))
								count += 1
								monomial_index[monomial] = count
								push!(model[:linear_operator], AffExpr(0., X_eq[l][i, j] => (2 - (i==j))*polynomial.coefficients[n])) 
							else
								add_to_expression!(model[:linear_operator][monomial_index[monomial]], (2 - (i==j))*polynomial.coefficients[n], X_eq[l][i, j])
							end

						end
						
					end
				end

			end

		end

		model[:X_eq] = X_eq

	end

	return nothing

end

function set_SOS_constraints!(model, pop, monomial_index)

	objective_coefficients = Dict{Int64, Float64}()

	for (support, coefficient) in terms(pop.objective)
		update!(objective_coefficients, monomial_index[support], coefficient)
	end

	add_to_expression!(model[:linear_operator][1], model[:t])

	set_variable_to_zero = AffExpr[]

	for (i, expression) in enumerate(model[:linear_operator])

		if i in keys(objective_coefficients)
			@constraint(model, expression == objective_coefficients[i])
		elseif length(expression.terms) == 1
			push!(set_variable_to_zero, first(expression.terms)[1])
		else
			@constraint(model, expression == 0.)
		end

	end

	for variable in unique(set_variable_to_zero)
		@constraint(model, variable == 0.)
	end

	return nothing

end

function sos_relaxation_model(pop::POP, relaxation_order::Int64, 
	variable_sets::Vector{Vector{T}},  
	# =dense_relaxation_set(pop),
	monomial_sets::Vector{Vector{Vector{UInt16}}},
	minimal_multipliers=false,
	return_monomials=false) where T <: Integer

	model = Model()

	@variable(model, t)

	monomial_index = set_X_0!(model, variable_sets, monomial_sets, relaxation_order)
	set_X_j!(model, pop, monomial_index, variable_sets, relaxation_order)
	set_SOS_constraints!(model, pop, monomial_index)

	@objective(model, Max, t)

	if return_monomials
		return model, monomial_index
	else
		return model
	end

end

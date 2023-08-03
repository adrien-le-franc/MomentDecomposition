# julia 1.6.5


function set_X_0!(model, sets, relaxation_order)

	monomial_index = Dict{Vector{UInt16}, Int64}()
	X_0 = Vector{Symmetric{VariableRef, Matrix{VariableRef}}}(undef, length(sets))
	linear_operator = Vector{AffExpr}()
	
	count = 0

	for (k, set) in enumerate(sets)

		matrix_size = n_moments(set, relaxation_order)
		X_0[k] = @variable(model, [1:matrix_size, 1:matrix_size], PSD)

		for (j, alpha) in moment_columns(set, relaxation_order)
			for (i, beta) in moment_rows(set, relaxation_order, j)

				monomial = monomial_product(alpha, beta)

				if !(monomial in keys(monomial_index))
					count += 1
					monomial_index[monomial] = count
					push!(linear_operator, AffExpr(0., X_0[k][i, j] => 2 - (i==j))) 
				else
					add_to_expression!(linear_operator[monomial_index[monomial]], 2 - (i==j), X_0[k][i, j])
				end

			end
		end

	end

	model[:X_0] = X_0
	model[:linear_operator] = linear_operator

	return monomial_index

end

function set_X_j!(model, pop, monomial_index, sets, relaxation_order)

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
					add_to_expression!(model[:linear_operator][monomial_index[support]], coefficient, X_ineq[l])
				end

			else

				matrix_size = n_moments(set, localizing_order)
				X_ineq[l] = @variable(model, [1:matrix_size, 1:matrix_size], PSD)

				for (j, alpha) in moment_columns(set, localizing_order)
					for (i, beta) in moment_rows(set, localizing_order, j)

						monomials = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support] 

						for (n, monomial) in enumerate(monomials)
							add_to_expression!(model[:linear_operator][monomial_index[monomial]], (2 - (i==j))*polynomial.coefficients[n], X_ineq[l][i, j])
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
					add_to_expression!(model[:linear_operator][monomial_index[support]], coefficient, X_eq[l])	
				end

			else

				matrix_size = n_moments(set, localizing_order)
				X_eq[l] = @variable(model, [1:matrix_size, 1:matrix_size], Symmetric)

				for (j, alpha) in moment_columns(set, localizing_order)
					for (i, beta) in moment_rows(set, localizing_order, j)

						monomials = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support] 

						for (n, monomial) in enumerate(monomials)
							add_to_expression!(model[:linear_operator][monomial_index[monomial]], (2 - (i==j))*polynomial.coefficients[n], X_eq[l][i, j])
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

	for (i, expression) in enumerate(model[:linear_operator])

		if i in keys(objective_coefficients)
			@constraint(model, expression == objective_coefficients[i])
		else
			@constraint(model, expression == 0.)
		end

	end

	return nothing

end

function sos_relaxation_model(pop::POP, relaxation_order::Int64, 
	sparsity_sets::Vector{Vector{T}}=dense_relaxation_set(pop), return_monomials=false) where T <: Integer

	model = Model()

	@variable(model, t)

	monomial_index = set_X_0!(model, sparsity_sets, relaxation_order)
	set_X_j!(model, pop, monomial_index, sparsity_sets, relaxation_order)
	set_SOS_constraints!(model, pop, monomial_index)

	@objective(model, Max, t)

	if return_monomials
		return model, monomial_index
	else
		return model
	end

end

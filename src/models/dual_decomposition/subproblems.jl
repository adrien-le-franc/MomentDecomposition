# julia 1.6
#
# dual decomposition subroblems for the hierarchy of moment relaxations


mutable struct SubProblem

	set_id::Int64
	model::Model

	# SubObjective ?
	f_k::Union{AffExpr, Nothing}
	coupling_terms::Union{Vector{Vector{AffExpr}}, Nothing} 

	# why ?
	multiplier_ids::Union{Vector{Int64}, Nothing}

	 
	
	#coupling_pairs::CouplingPairs 

	#moment_labels::Union{Dict{Vector{UInt16}, Int64}, Nothing} # besoin ?

end

function Subproblem(k::Int64)
	return Subproblem(k, Model(), nothing, nothing, nothing)
end

### ?
struct MultiplierInformation
	moments::Vector{Vector{UInt16}}
	pair::Tuple{UInt64, UInt16}
end




function set_polynomial_objective!(subproblems::Vector{Subproblem}, pop::POP, 
	relaxation_order::Int64, moment_labels, variable_sets::Vector{Vector{UInt16}})
	
	objective = objective_decomposition(pop, variable_sets) ## allocations -> slow ??

	for k in 1:length(variable_sets)

		n_temrs = length(objective[k].coefficients)
		
		if n_temrs > 0

			moments_f_k = [moment_labels[k][objective[k].support[n]] for n in 1:n_temrs]
			
			subproblems[k].f_k = @expression(subproblems[k].model, 
				sum(objective[k].coefficients[n]*subproblems[k].model[Symbol("y_$k")][moments_f_k[n]] for n in 1:n_temrs))
		
		end
	

	end

	return nothing

end

function set_inequality_constraints!(subproblems::Vector{Subproblem}, pop::POP, 
	relaxation_order::Int64, moment_labels, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.inequality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(subproblems[k].model, 
				sum(polynomial.coefficients[n]*subproblems[k].models[Symbol("y_$k")][moment_labels[k][polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) >= 0)

		else

			n_monomials = n_moments(length(variable_sets[k]), relaxation_order)
			localizing_matrix = @variable(subproblems[k].model, [1:n_monomials, 1:n_monomials], PSD)

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(subproblems[k].model, 
						localizing_matrix[i, j] == sum(polynomial.coefficients[n]*subproblems[k].models[Symbol("y_$k")][moment_labels[k][labels[n]]] for n in 1:length(polynomial.coefficients)))

				end
			end

		end


	end

	return nothing

end

function set_equality_constraints!(subproblems::Vector{Subproblem}, pop::POP, 
	relaxation_order::Int64, moment_labels, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.equality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(subproblems[k].model, 
				sum(polynomial.coefficients[n]*subproblems[k].model[Symbol("y_$k")][moment_labels[k][polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) == 0)

		else

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(subproblems[k].model, 
						sum(polynomial.coefficients[n]*subproblems[k].model[Symbol("y_$k")][moment_labels[k][labels[n]]] for n in 1:length(polynomial.coefficients)) == 0.)

				end
			end			

		end


	end

	return nothing

end

function set_Lagrange_multipliers!(subproblems::Vector{Subproblem}, 
	relaxation_order, moment_labels, variable_sets, max_coupling_order)



	multiplier_information = MultiplierInformation[] 


	for (id, pair) in enumerate(combinations(1:length(variable_sets), 2))

		k_1, k_2 = pair
		intersection = intersect(variable_sets[k_1], variable_sets[k_2])

		if isempty(intersection)
			continue
		end

		pair_information = MultiplierInformation(Vector{UInt16}[] , (k_1, k_2))

		coupling_moments_k_1 = AffExpr[]
		coupling_moments_k_2 = AffExpr[]

		for alpha in moments(intersection, 2*relaxation_order) # remove moment 0...0 ?

			label = monomial(alpha)

			if length(label) > max_coupling_order
				continue
			end

			# store coupling information in subproblems

			push!(coupling_moments_k_1, 1*subproblems[k_1].model[Symbol("y_$(k_1)")][moment_labels[k_1][label]])
			push!(coupling_moments_k_2, -1*subproblems[k_2].model[Symbol("y_$(k_2)")][moment_labels[k_2][label]])

			# store coupling information in master

			push!(pair_information.moments, label)

		end

		push!(subproblems[k_1].coupling_terms, coupling_moments_k_1)
		push!(subproblems[k_2].coupling_terms, coupling_moments_k_2)
		push!(subproblems[k_1].multiplier_ids, id)
		push!(subproblems[k_2].multiplier_ids, id)

		push!(multiplier_information, pair_information)

	end

	return multiplier_information

end

function update_dual_objective!(subproblem::SubProblem, submultiplier::Vector{Vector{Float64}})

	# normalize ?

	@objective(subproblem.model, Min, 
		subproblem.f_k + sum(subproblem.coupling_terms[k]'*submultiplier[k] 
		for k in 1:length(submultiplier)))

	return nothing

end

function subproblem_oracle(subproblem, submultiplier)

	update_dual_objective!(subproblem, submultiplier)
	optimize!(subproblem.model)


	# check optimizer status !!

	return push!([value.(coupling_term) for coupling_term in subproblem.coupling_terms], 
			[objective_value(subproblem.model)]) # check super/sub gradient !!

end
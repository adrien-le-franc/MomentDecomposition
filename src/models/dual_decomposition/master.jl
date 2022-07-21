# julia 1.6
#
# dual decomposition master problem for the hierarchy of moment relaxations


function dual_decomposed_relaxation(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
	max_coupling_order::Int64=maxInt64)

	subproblems, coupling_moments = build_dual_decomposition(pop, relaxation_order, 
											variable_sets, max_coupling_order)



	optimize_dual_decomposition!()



end

function build_dual_decomposition(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
	max_coupling_order::Int64=maxInt64) 

	#n_sets = length(variable_sets)
	
	# distribute f_k
	objective = objective_decomposition(pop, variable_sets)


	n_sets = length(variable_sets)
	subproblems = [Subproblem(k, objective[k]) for k in 1:n_sets]

	# build subproblems : moment variables and labels, moment matrices, probability constraint
	## run in parallel ?

	for (k, set) in enumerate(variable_sets)

		set_moment_variables!(subproblems[k], relaxation_order, set, k)
		set_moment_matrix!(subproblems[k], relaxation_order, set, k)
		set_probability_measure_constraint!(subproblems[k], set, k) 

	end

	# distribute polynomials over variables sets
	## run in parallel ?

	#set_objective!(model, pop, moment_labels, variable_sets)

	set_inequality_constraints!(subproblems, pop, relaxation_order, variable_sets)
	set_equality_constraints!(subproblems, pop, relaxation_order, variable_sets)


	# set dual variable / coupling variable

	coupling_moments = set_coupling_constraints!(model, relaxation_order, 
							moment_labels, variable_sets, max_coupling_order)

	
	return subproblems, coupling_moments

end
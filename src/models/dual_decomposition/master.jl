# julia 1.6
#
# dual decomposition master problem for the hierarchy of moment relaxations


function build_dual_subproblems(variable_sets::Vector{Vector{T}}, 
	relaxation_order::Int64) where T<:Integer

	n_sets = length(variable_sets)
	subproblems = [Subproblem(k) for k in 1:n_sets]
	moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:n_sets)

	for (k, set) in enumerate(variable_sets)

		set_moment_variables!(subproblems[k].model, set, relaxation_order)
		moment_labels[k] = set_moment_matrix!(subproblems[k].model, set, relaxation_order)
		set_probability_measure_constraint!(subproblems[k].model, moment_labels[k]) 

	end

	return subproblems, moment_labels

end

function dual_decomposition_models(pop::POP, relaxation_order::Int64, 
	variable_sets::Vector{Vector{T}}; max_coupling_order::Int64=maxInt64) where T<:Integer 
	
	subproblems, moment_labels = build_dual_subproblems(variable_sets, relaxation_order)

	set_polynomial_objective!(subproblems, pop, moment_labels, variable_sets)
	set_polynomial_constraints!(subproblems, pop, relaxation_order, moment_labels, variable_sets)
	
	multiplier_information = set_Lagrange_multipliers!(subproblems, relaxation_order, 
								moment_labels, variable_sets, max_coupling_order)

	return subproblems, multiplier_information

end
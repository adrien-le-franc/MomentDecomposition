# julia 1.6
#
# dual decomposition master problem for the hierarchy of moment relaxations


mutable struct Multiplier
	value::Vector{Vector{Float64}}
	information::Vector{MultiplierInformation}
	#n_pairs::Int64
end



function build_dual_decomposition(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
	max_coupling_order::Int64=maxInt64) 
	
	# init moment labels
	n_sets = length(variable_sets)
	moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:n_sets)

	# init subproblems
	subproblems = [Subproblem(k) for k in 1:n_sets]

	# build subproblems : moment variables and labels, moment matrices, probability constraint
	## run in parallel ?

	for (k, set) in enumerate(variable_sets)

		set_moment_variables!(subproblems[k].model, relaxation_order, set, k)
		moment_labels[k] = set_moment_matrix!(subproblems[k].model, relaxation_order, set, k)
		set_probability_measure_constraint!(subproblems[k].model, moment_labels, set, k) 

	end

	# distribute polynomials over variables sets
	## run in parallel ?

	set_polynomial_objective!(subproblems, pop, relaxation_order, moment_labels, variable_sets)
	set_inequality_constraints!(subproblems, pop, relaxation_order, moment_labels, variable_sets)
	set_equality_constraints!(subproblems, pop, relaxation_order, moment_labels, variable_sets)

	# set lagrange multipliers

	multiplier_information = set_Lagrange_multipliers!(subproblems, relaxation_order, 
							moment_labels, variable_sets, max_coupling_order)

	
	return subproblems, multiplier_information

end

function call_maximization_oracle!(subproblems, multiplier, oracle_data)

	oracle_data = pmap(subproblem->subproblem_oracle(subproblem, 
		extract(multiplier, subproblem)), subproblems)

	return supergradient(oracle_data, multiplier), dual_objective_value(oracle_data)

end

function extract(multiplier::Multiplier, subproblem::SubProblem)

	return [multiplier.value[i] for i in subproblem.multiplier_ids] ### order ??????? ou alors sort ?

end

function supergradient(subproblems, multiplier, oracle_data)

	return

	[
		oracle_data[information.pair[1]][nbcfind(subproblems[information.pair[1]].multiplier_ids, k)] +
		oracle_data[information.pair[2]][nbcfind(subproblems[information.pair[2]].multiplier_ids, k)] 

		for (k, information) in enumerate(multiplier.information)
	]

end

function dual_objective_value(oracle_data)
	return sum(subproblem_data[end][1] for subproblem_data in oracle_data)
end

function initialize_multiplier(multiplier_value, multiplier_information)

	return ?

end


function run_optimization!(subproblems::Vector{SubProblem}, multiplier::Multiplier)

	# ??? update Subgradient Methods first

end



function dual_decomposed_relaxation(pop::POP, relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}};
	max_coupling_order::Int64=maxInt64, 
	multiplier_value::Vector{Vector{Float64}}=[zeros(length(set)) for set in variable_sets])

	subproblems, multiplier_information = build_dual_decomposition(pop, relaxation_order, 
											variable_sets, max_coupling_order)

	multiplier = Multiplier(multiplier_value, multiplier_information)

 	run_optimization!(subproblems, multiplier)

end
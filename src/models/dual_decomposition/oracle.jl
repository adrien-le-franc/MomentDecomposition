# julia 1.6
#
# dual decomposition oracle for the hierarchy of moment relaxations

function update_dual_objective!(subproblem::SubProblem, submultiplier::Vector{Vector{Float64}})

	# normalize ?

	@objective(subproblem.model, Min, 
		subproblem.polynomial_objective + sum(subproblem.coupling_terms[k]'*submultiplier[k] 
		for k in 1:length(submultiplier)))

	return nothing

end

function call_subproblem_oracle!(subproblem::SubProblem, submultiplier::Vector{Vector{Float64}})

	update_dual_objective!(subproblem, submultiplier)
	optimize!(subproblem.model)

	# check optimizer status !!

	return push!([value.(coupling_term) for coupling_term in subproblem.coupling_terms], [objective_value(subproblem.model)]) # check super/sub gradient !!

end

function extract(multiplier::Multiplier, subproblem::SubProblem)

	return [multiplier.value[i] for i in subproblem.multiplier_ids]

end

function supergradient(subproblems::Vector{SubProblem}, multiplier::Multiplier, 
	oracle_data::Vector{Vector{Vector{Float64}}})

	return [
		oracle_data[information.pair[1]][searchsortedfirst(subproblems[information.pair[1]].multiplier_ids, k)] +
		oracle_data[information.pair[2]][searchsortedfirst(subproblems[information.pair[2]].multiplier_ids, k)] 

		for (k, information) in enumerate(multiplier.information)]

end

function dual_objective_value(oracle_data::Vector{Vector{Vector{Float64}}})
	return sum(subproblem_data[end][1] for subproblem_data in oracle_data)
end

function call_maximization_oracle!(subproblems::Vector{SubProblem}, 
	multiplier::Multiplier, oracle_data::Vector{Vector{Vector{Float64}}})

	oracle_data = pmap(subproblem->call_subproblem_oracle!(subproblem, 
		extract(multiplier, subproblem)), subproblems)

	return nothing 

end

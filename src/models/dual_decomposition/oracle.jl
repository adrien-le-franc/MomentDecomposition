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

function call_subproblem_oracle!(subproblem, submultiplier)

	update_dual_objective!(subproblem, submultiplier)
	optimize!(subproblem.model)


	# check optimizer status !!

	return vcat([value.(coupling_term) for coupling_term in subproblem.coupling_terms], 
			[objective_value(subproblem.model)]) # check super/sub gradient !!

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

function call_maximization_oracle!(subproblems, multiplier, oracle_data)

	oracle_data = pmap(subproblem->subproblem_oracle(subproblem, 
		extract(multiplier, subproblem)), subproblems)

	return supergradient(oracle_data, multiplier), dual_objective_value(oracle_data)

end

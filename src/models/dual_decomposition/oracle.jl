# julia 1.6
#
# dual decomposition oracle for the hierarchy of moment relaxations

function update_dual_objective!(subproblem::SubProblem, submultiplier::Vector{Vector{Float64}})

	scale_factor = maximum(vcat(abs.(values(subproblem.polynomial_objective.terms)), 
								[maximum(abs.(l)) for l in submultiplier]))

	if scale_factor == 0.
		@objective(subproblem.model, Min, 0.)
	else
		@objective(subproblem.model, Min, 
			(subproblem.polynomial_objective + sum(subproblem.coupling_terms[k]'*submultiplier[k] 
			for k in 1:length(submultiplier))) / scale_factor )
	end

	return scale_factor

end

function call_subproblem_oracle!(subproblem::SubProblem, submultiplier::Vector{Vector{Float64}})

	scale_factor = update_dual_objective!(subproblem, submultiplier)
	optimize!(subproblem.model)

	# check optimizer status !!
	"""
	println(solution_summary(subproblem.model))
	println(termination_status(subproblem.model))
	println(primal_status(subproblem.model))
	"""

	if primal_status(subproblem.model) != FEASIBLE_POINT
		println("WARNING: solution returned is primal unfeasible")
	end

	return push!([value.(coupling_term) for coupling_term in subproblem.coupling_terms], 
		[objective_value(subproblem.model) * scale_factor]) # check super/sub gradient !!

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

	oracle_data[:] = pmap(subproblem->call_subproblem_oracle!(subproblem, # perf ???
		extract(multiplier, subproblem)), subproblems)

	return nothing 

end

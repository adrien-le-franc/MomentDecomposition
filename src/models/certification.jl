# julia 1.6
#
# certified bounds for SOS relaxations


function project_to_PSD_cone(M::Symmetric{Float64, Matrix{Float64}})

	sigmas, P = eigen(M, sortby=x->-x) # allocations !? is sorting necessary ?
	return Symmetric(P*spdiagm(max.(sigmas, 0.))*P')

end

function project_to_PSD_cone(X::Symmetric{VariableRef, Matrix{VariableRef}})
	return project_to_PSD_cone(Symmetric(value.(X)))
end

function project_to_PSD_cone(x::VariableRef)
	return max(0., value(x))
end

function get_solver_values(model::Model, pop::POP)

	values = Dict{VariableRef, Float64}()

	for k in 1:length(model[:X_0])
		X = project_to_PSD_cone(model[:X_0][k])
		n = size(X)[1]
		for i in 1:n
			for j in 1:i
				values[model[:X_0][k][i, j]] = X[i, j]
			end
		end
	end

	for l in 1:length(pop.inequality_constraints)
		X = project_to_PSD_cone(model[:X_ineq][l])
		if typeof(X) == Float64
			values[model[:X_ineq][l]] = X
		else
			n = size(X)[1]
			for i in 1:n
				for j in 1:i
					values[model[:X_ineq][l][i, j]] = X[i, j]
				end
			end	
		end
	end

	for l in 1:length(pop.equality_constraints)
		X = value.(model[:X_eq][l])
		if typeof(X) == Float64
			values[model[:X_eq][l]] = X
		else
			n = size(X)[1]
			for i in 1:n
				for j in 1:i
					values[model[:X_eq][l][i, j]] = X[i, j]
				end
			end	
		end
	end

	values[model[:t]] = value(model[:t])

	return values

end

function update!(dict::Dict{Int64, Float64}, index::Int64, coefficient::Float64)

	if index in keys(dict)
		dict[index] += coefficient
	else
		dict[index] = coefficient
	end

	return nothing
end

function evaluate_linear_system_error(model::Model, pop::POP, monomial_index::Dict{Vector{UInt16}, Int64},
	values::Dict{VariableRef, Float64})

	objective_coefficients = Dict{Int64, Float64}()

	for (support, coefficient) in terms(pop.objective)
		update!(objective_coefficients, monomial_index[support], coefficient)
	end

	n_monomials = length(model[:linear_operator])
	err = zeros(n_monomials)

	for (i, expression) in enumerate(model[:linear_operator])

		if i in keys(objective_coefficients)
			err[i] += value(z->values[z], expression) - objective_coefficients[i]
		else
			err[i] += value(z->values[z], expression)
		end

	end

	return sum([e <= 0. ? e : 0. for e in -1*err])

end

function compute_error_for_scaled_sos(model::Model, pop::POP, monomial_index::Dict{Vector{UInt16}, Int64})

	values = get_solver_values(model, pop)
	return evaluate_linear_system_error(model, pop, monomial_index, values)

end
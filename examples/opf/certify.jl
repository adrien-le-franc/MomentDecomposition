# julia 1.6.5
#
# bound certifyer

using LinearAlgebra
using SparseArrays

function projection(M::Symmetric{Float64, Matrix{Float64}}) # keep it sparse untill here ? Symmetric after Matrix ?

	sigmas, P = eigen(M, sortby=x->-x) # allocations !? is sorting necessary ?

	return Symmetric(P*spdiagm(max.(sigmas, 0.))*P')

end

function proj(X::Symmetric{VariableRef, Matrix{VariableRef}})
	return projection(Symmetric(value.(X)))
end

function proj(x::VariableRef)
	return max(0., value(x))
end

function compute_error(model, pop, monomial_index)

	## project X

	values = Dict{VariableRef, Float64}()

	# sigma 0
	for k in 1:length(model[:X_0])
		X = proj(model[:X_0][k])
		n = size(X)[1]
		for i in 1:n
			for j in 1:i
				values[model[:X_0][k][i, j]] = X[i, j]
			end
		end
	end

	# sigma ineq
	for l in 1:length(pop.inequality_constraints)
		X = proj(model[:X_ineq][l])
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

	# sigma eq
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

	# t
	values[model[:t]] = value(model[:t])

	## eval error

	objective_coefficients = Dict{Int64, Float64}()

	for (support, coefficient) in MSOS.terms(pop.objective)
		update!(objective_coefficients, monomial_index[support], coefficient)
	end

	n_monomials = length(monomial_index)
	err = zeros(n_monomials)

	for (i, expression) in enumerate(model[:linear_operator])

		if i in keys(objective_coefficients)
			err[i] += value(z->values[z], expression) - objective_coefficients[i]
		else
			err[i] += value(z->values[z], expression)
		end

	end

	return -1*err

end
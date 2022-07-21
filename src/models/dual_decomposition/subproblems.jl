# julia 1.6
#
# dual decomposition subroblems for the hierarchy of moment relaxations


mutable struct Subproblem

	set_id::Int64
	model::Model

	coupling_terms::Vector{Vector{AffExpr}} #### best deal ! Dict ?

	f_k::SparsePolynomial ### AffExpr
	
	#coupling_pairs::CouplingPairs 

	#moment_labels::Union{Dict{Vector{UInt16}, Int64}, Nothing} # besoin ?

end

function Subproblem(k::Int64, f::SparsePolynomial)
	return Subproblem(k, Model(), f, CouplingPairs(), nothing)
end

struct CouplingPairs
	ranks::Vector{Int64}
	signs::Vector{Float64}
end

function CouplingPairs() = CouplingPairs(Int64[], Float64[])

struct CouplingMoments
	moments::Vector{Vector{UInt16}}
	pair::Tuple{UInt64, UInt16}
end

function set_moment_variables!(subproblem::Subproblem, relaxation_order::Int64, 
	set::Vector{UInt16}, k::Int64)

	set_moment_variables!(subproblem.model, relaxation_order, set, k)

	return nothing

end

function set_moment_matrix!(subproblem::Subproblem, relaxation_order::Int64, set::Vector{UInt16}, k::Int64)
	
	subproblem.moment_labels = set_moment_matrix!(subproblem.model, relaxation_order, set, k)

	return nothing

end

function set_probability_measure_constraint!(subproblem::Subproblem, set, k)
	@constraint(subproblem.model, subproblem.model[Symbol("y_$k")][subproblem.moment_labels[UInt16[]]] == 1.)
	return nothing
end


function set_inequality_constraints!(subproblems::Vector{Subproblem}, pop::POP, 
	relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.inequality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(subproblems[k].model, 
				sum(polynomial.coefficients[n]*subproblems[k].models[Symbol("y_$k")][subproblems[k].moment_labels[polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) >= 0)

		else

			n_monomials = n_moments(length(variable_sets[k]), relaxation_order)
			localizing_matrix = @variable(subproblems[k].model, [1:n_monomials, 1:n_monomials], PSD)

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(subproblems[k].model, 
						localizing_matrix[i, j] == sum(polynomial.coefficients[n]*subproblems[k].models[Symbol("y_$k")][subproblems[k].moment_labels[labels[n]]] for n in 1:length(polynomial.coefficients)))

				end
			end

		end


	end

	return nothing

end

function set_equality_constraints!(subproblems::Vector{Subproblem}, pop::POP, 
	relaxation_order::Int64, variable_sets::Vector{Vector{UInt16}})

	for polynomial in pop.equality_constraints

		k = assign_constraint_to_set(polynomial, variable_sets)
		localizing_order = localizing_matrix_order(relaxation_order, polynomial)

		if localizing_order == 0

			@constraint(subproblems[k].model, 
				sum(polynomial.coefficients[n]*subproblems[k].model[Symbol("y_$k")][subproblems[k].moment_labels[polynomial.support[n]]]
				for n in 1:length(polynomial.coefficients)) == 0)

		else

			for (j, alpha) in moment_columns(variable_sets[k], localizing_order)
				for (i, beta) in moment_rows(variable_sets[k], localizing_order, j)

					labels = [monomial_product(alpha, beta, gamma) for gamma in polynomial.support]
					@constraint(subproblems[k].model, 
						sum(polynomial.coefficients[n]*subproblems[k].model[Symbol("y_$k")][subproblems[k].moment_labels[labels[n]]] for n in 1:length(polynomial.coefficients)) == 0.)

				end
			end			

		end


	end

	return nothing

end

function set_coupling_constraints!(subproblems::Vector{Subproblem}, 
	relaxation_order, moment_labels, variable_sets, max_coupling_order)

	all_coupling_moments = CouplingMoments[]

	for (i, pair) in enumerate(combinations(1:length(variable_sets), 2))

		k_1, k_2 = pair
		intersection = intersect(variable_sets[k_1], variable_sets[k_2])

		if isempty(intersection)
			continue
		end

		coupling_moments = CouplingMoments(Vector{UInt16}[] , (k_1, k_2))

		for alpha in moments(intersection, 2*relaxation_order) # remove moment 0...0 ?

			label = monomial(alpha)

			if length(label) > max_coupling_order
				continue
			end

			# store coupling information in subproblems

			push!(subproblems[k_1].coupling_pairs.ranks, i)
			push!(subproblems[k_1].coupling_pairs.signs, 1.)

			push!(subproblems[k_2].coupling_pairs.ranks, i)
			push!(subproblems[k_2].coupling_pairs.signs, -1.)

			# store coupling information in master

			push!(coupling_moments.moments, label)

		end

		push!(all_coupling_moments, coupling_moments)

	end

	return all_coupling_moments

end

function set_dual_objective!(subproblem::Subproblem, moment_labels, dual_variable)


	f = AffExpr(0.) 

	# objective f_k -> def as an expression once and for all ???

	n_temrs = length(subproblem.f_k.coefficients)
	
	if n_temrs > 0
		
		moments_f_k = [subproblem.moment_labels[subproblem.f_k.support[n]] for n in 1:n_temrs]
		f += sum(subproblem.f_k.coefficients[n]*subproblem.model[Symbol("y_$k")][moments_f_k[n]] for n in 1:n_temrs)	
	
	end

	# coupling terms -> def moment vector as an expression once and for all ???






	@objective(model, Min, f)

	return nothing

end
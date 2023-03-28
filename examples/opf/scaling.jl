# julia 1.6

using MomentSOS
MSOS = MomentSOS

using Combinatorics


function update!(epsilon, support, coefficient)

	if support in keys(epsilon)
		epsilon[support] += coefficient
	else
		epsilon[support] = coefficient
	end

	return nothing
end

function dict_to_polynomial(dict::Dict{Vector{UInt16}, Float64})

	support = Vector{UInt16}[]
	coefficient = Float64[]

	for (key, val) in dict

		push!(support, key)
		push!(coefficient, val)

	end

	return MSOS.SparsePolynomial(support, coefficient)

end

function scale_polynomial(polynomial::MSOS.SparsePolynomial, bounds::Dict{UInt16, Vector{Float64}}, normalize::Bool=false)

	p = Dict{Vector{UInt16}, Float64}()

	for (support, coefficient) in MSOS.terms(polynomial)

		variables = unique(support)

		if length(variables) == 0
			update!(p, support, coefficient)
		else

			degrees = [sum(support .== v) for v in variables]

			for alpha in Iterators.product([0:d for d in degrees]...)

				new_support = UInt16[]
				new_coefficient = coefficient

				for (i, k_i) in enumerate(alpha)

					append!(new_support, [variables[i] for _ in 1:k_i])

					u = bounds[variables[i]][1] - bounds[variables[i]][2]
					v = bounds[variables[i]][2]
					new_coefficient = new_coefficient*u^k_i*v^(degrees[i] - k_i)*binomial(degrees[i], k_i)

				end

				update!(p, new_support, new_coefficient)

			end

		end

	end

	if normalize

		max_coefficient = maximum(abs.(values(p)))
		support = Vector{UInt16}[]
		coefficient = Float64[]

		for (key, val) in p

			if abs(val / max_coefficient) >= 1e-8

				push!(support, key)
				push!(coefficient, val / max_coefficient)

			end

		end

		return max_coefficient, MSOS.SparsePolynomial(support, coefficient)
	
	else
		
		return dict_to_polynomial(p)
	
	end

end

function normalize_polynomial(p::MSOS.SparsePolynomial)

	max_coefficient = maximum(abs.(p.coefficients))
	
	support = Vector{UInt16}[]
	coefficient = Float64[]
	new_c = 0.

	for (s, c) in MSOS.terms(p)

		new_c = c / max_coefficient

		if abs(new_c) >= 1e-8
			push!(support, s)
			push!(coefficient, new_c)
		end

	end

	return MSOS.SparsePolynomial(support, coefficient)

end


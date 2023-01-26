# julia 1.6.5

# same types for equalities and inequalities may be a bad idea...

mutable struct Phi_scalar
	f::Function
	f_and_g!::Function
end

mutable struct Phi_matrix 
	M::Symmetric{Float64, Matrix{Float64}} # change for sparse matrix in equality constraints ?
	A_j::Vector{SparseMatrixCSC{Float64, Int64}}
	address::Vector{Int64}
end

mutable struct Phi_localizing
	scalar::Vector{Phi_scalar}
	matrix::Vector{Phi_matrix}
end

mutable struct DualModel
	phi_moment::Phi_matrix
	phi_inequality::Phi_localizing
	phi_equality::Phi_localizing
	phi_0::Phi_scalar
	linear_term::Phi_scalar
end


function set_phi_moment(pop::POP, relaxation_order::Int64)

	moment_labels = Dict{Vector{UInt16}, Int64}()
	matrix_size = n_moments(pop, relaxation_order)
	entries_A_j = Dict{UInt16, Matrix{Float64}}()

	count = 0

	for (j, alpha) in moment_columns(pop, relaxation_order)
		for (i, beta) in moment_rows(pop, relaxation_order, j)

			label = monomial_product(alpha, beta)

			if !(label in keys(moment_labels))
				count += 1
				moment_labels[label] = count
				entries_A_j[count] = [matrix_size matrix_size 0.; i j 1.]
			else
				entries_A_j[moment_labels[label]] = vcat(entries_A_j[moment_labels[label]], [i j 1.])
			end

		end
	end

	A_j = [dropzeros(sparse(convert.(Int64, entries_A_j[i][:, 1]), 
				  			convert.(Int64, entries_A_j[i][:, 2]),
			  	  							entries_A_j[i][:, 3])) for i in 1:count]

	M = Symmetric(Matrix{Float64}(undef, matrix_size, matrix_size))

	return Phi_matrix(M, A_j, collect(1:count)), moment_labels

end

function set_phi_matrix(polynomial::SparsePolynomial, pop::POP, localizing_order::Int64, 
	moment_labels::Dict{Vector{UInt16}, Int64})

	n_monomials = n_moments(pop, localizing_order)
	entries_A_j = Dict{Vector{UInt16}, Matrix{Float64}}()

	matrix_size = n_moments(pop, localizing_order)
	entries_A_j = Dict{UInt16, Matrix{Float64}}()

	for (j, alpha) in moment_columns(pop, localizing_order)
		for (i, beta) in moment_rows(pop, localizing_order, j)

			labels = [moment_labels[monomial_product(alpha, beta, gamma)] for gamma in polynomial.support] # confusion : labels ??? alpha or i ?

			for (k, label) in enumerate(labels)

				if !(label in keys(entries_A_j))
					entries_A_j[label] = [matrix_size matrix_size 0.; i j polynomial.coefficients[k]]
				else
					entries_A_j[label] = vcat(entries_A_j[label], [i j polynomial.coefficients[k]])
				end	

			end
			
		end
	end

	address = sort!(collect(keys(entries_A_j)))
	A_j = [dropzeros(sparse(convert.(Int64, entries_A_j[i][:, 1]), 
				  			convert.(Int64, entries_A_j[i][:, 2]),
			  	  							entries_A_j[i][:, 3])) for i in address] 

	M = Symmetric(Matrix{Float64}(undef, n_monomials, n_monomials))

	return Phi_matrix(M, A_j, address)

end

function get_linear_form(polynomial::SparsePolynomial, moment_labels::Dict{Vector{UInt16}, Int64})

	entries_a_j = Dict{UInt16, Float64}()

	for (support, coefficient) in terms(polynomial)

		label = moment_labels[support]

		if !(label in keys(entries_a_j))
			entries_a_j[label] = coefficient
		else
			entries_a_j[label] += coefficient
		end

	end

	address = sort!(collect(keys(entries_a_j))) 

	return [entries_a_j[i] for i in address], address

end

function set_phi_scalar_inequality(polynomial::SparsePolynomial, pop::POP, moment_labels::Dict{Vector{UInt16}, Int64})

	a_j, address = get_linear_form(polynomial, moment_labels)

	@views f(y::Vector{Float64}) = -0.5*max(a_j'*y[address], 0.)^2 

	@views function f_and_g!(y::Vector{Float64}, grad::Vector{Float64}) 

		max_ay = max(a_j'*y[address], 0)
		grad[address] += -a_j*max_ay

		return -0.5*max_ay^2

	end

	return Phi_scalar(f, f_and_g!)

end

function set_phi_scalar_equality(polynomial::SparsePolynomial, pop::POP, moment_labels::Dict{Vector{UInt16}, Int64})

	a_j, address = get_linear_form(polynomial, moment_labels)

	@views f(y::Vector{Float64}) = -0.5*(a_j'*y[address])^2 

	@views function f_and_g!(y::Vector{Float64}, grad::Vector{Float64}) 

		ay = a_j'*y[address]
		grad[address] += -a_j*ay

		return -0.5*ay^2

	end

	return Phi_scalar(f, f_and_g!)

end

function set_phi_inequality(pop::POP, relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64})

	phi_matrix = Vector{Phi_matrix}()
	phi_scalar = Vector{Phi_scalar}()

	if pop.inequality_constraints == nothing
		return Phi_localizing(phi_scalar, phi_matrix)
	else
		for polynomial in pop.inequality_constraints

			localizing_order = localizing_matrix_order(polynomial, relaxation_order)

			if localizing_order == 0
				push!(phi_scalar, set_phi_scalar_inequality(polynomial, pop, moment_labels))
			else
				push!(phi_matrix, set_phi_matrix(polynomial, pop, localizing_order, moment_labels)) 
			end

		end
	end

	return Phi_localizing(phi_scalar, phi_matrix)

end

function set_phi_equality(pop::POP, relaxation_order::Int64, moment_labels::Dict{Vector{UInt16}, Int64})

	phi_matrix = Vector{Phi_matrix}()
	phi_scalar = Vector{Phi_scalar}()

	if pop.equality_constraints == nothing
		return Phi_localizing(phi_scalar, phi_matrix)
	else
		for polynomial in pop.equality_constraints

			localizing_order = localizing_matrix_order(polynomial, relaxation_order)

			if localizing_order == 0
				push!(phi_scalar, set_phi_scalar_equality(polynomial, pop, moment_labels))
			else
				push!(phi_matrix, set_phi_matrix(polynomial, pop, localizing_order, moment_labels)) 
			end

		end
	end

	return Phi_localizing(phi_scalar, phi_matrix)

end

function set_phi_0(value_to_certify::Float64)

	@views function f(y::Vector{Float64})

		if y[1] > 0
			return -0.5*y[1]^2 - value_to_certify*y[1]
		else
			return -value_to_certify*y[1]
		end

	end

	@views function f_and_g!(y::Vector{Float64}, grad::Vector{Float64})

		if y[1] > 0
			grad[1] += -y[1] - value_to_certify
			return -0.5*y[1]^2 - value_to_certify*y[1] 
		else
			grad[1] += -value_to_certify
			return -value_to_certify*y[1] 
		end

	end

	return Phi_scalar(f, f_and_g!)

end

function set_linear_term(pop::POP, moment_labels::Dict{Vector{UInt16}, Int64})

	a_j, address = get_linear_form(pop.objective, moment_labels)

	@views f(y::Vector{Float64}) = a_j'*y[address]

	@views function f_and_g!(y::Vector{Float64}, grad::Vector{Float64})

		grad[address] += a_j

		return a_j'*y[address]

	end

	return Phi_scalar(f, f_and_g!)

end

function build_dual_model(pop::POP, relaxation_order::Int64, value_to_certify::Float64)

	phi_moment, moment_labels = set_phi_moment(pop, relaxation_order)
	phi_inequality = set_phi_inequality(pop, relaxation_order, moment_labels)
	phi_equality = set_phi_equality(pop, relaxation_order, moment_labels)

	phi_0 = set_phi_0(value_to_certify)
	phi_linear_term = set_linear_term(pop, moment_labels)

	return DualModel(phi_moment, phi_inequality, phi_equality, phi_0, phi_linear_term)

end


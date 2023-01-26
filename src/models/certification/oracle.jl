# julia 1.6.5

## oracle ## -> improve solver flexibility (common computations, f, g, fg ...)

#call_f!(y::Vector{Float64}, phi::Phi_scalar) = phi.f(y)

call_f(y::Vector{Float64}, phi::Phi_scalar) = phi.f(y)
call_f_and_g!(y::Vector{Float64}, grad::Vector{Float64}, phi::Phi_scalar) = phi.f_and_g!(y, grad)

@views function evaluate_matrix(y::Vector{Float64}, phi::Phi_matrix) 
	return Symmetric(Matrix(sum(phi.A_j[i]*y[phi.address][i] for i in 1:length(phi.A_j)))) # keep sparse matrix for equality constraints ?
end

function projection(M::Symmetric{Float64, Matrix{Float64}}) # keep it sparse untill here ? Symmetric after Matrix ?

	sigmas, P = eigen(M, sortby=x->-x) # allocations !? is sorting necessary ?

	return Symmetric(P*spdiagm(max.(sigmas, 0.))*P')

end

function update_matrices!(y::Vector{Float64}, dual_model::DualModel)

	dual_model.phi_moment.M = projection(evaluate_matrix(y, dual_model.phi_moment))

	for phi in dual_model.phi_inequality.matrix
		phi.M = projection(evaluate_matrix(y, phi))
	end

	for phi in dual_model.phi_equality.matrix
		phi.M = evaluate_matrix(y, phi)
	end

	return nothing

end

function call_g(phi::Phi_matrix)
	return -1*[dot(Symmetric(phi.A_j[i]), phi.M) for i in 1:length(phi.A_j)]
end

function update_gradient!(gradient::Vector{Float64}, y::Vector{Float64}, dual_model::DualModel)

	gradient[:] = call_g(dual_model.phi_moment)

	for phi in dual_model.phi_inequality.matrix
		gradient[phi.address] += call_g(phi)
	end

	for phi in dual_model.phi_inequality.scalar
		_ = call_f_and_g!(y, gradient, phi) # update syntax ?
	end

	for phi in dual_model.phi_equality.matrix
		gradient[phi.address] += call_g(phi)
	end

	for phi in dual_model.phi_equality.scalar
		_ = call_f_and_g!(y, gradient, phi)
	end	

	_ = call_f_and_g!(y, gradient, dual_model.phi_0)
	_ = call_f_and_g!(y, gradient, dual_model.linear_term)

	return nothing

end

function call_f(phi::Phi_matrix)
	return -0.5*norm(phi.M)^2
end

function evaluate_objective(y::Vector{Float64}, dual_model::DualModel)

	obj = call_f(dual_model.phi_moment)

	for phi in dual_model.phi_inequality.matrix
		obj += call_f(phi)
	end

	for phi in dual_model.phi_inequality.scalar
		obj += call_f(y, phi)
	end

	for phi in dual_model.phi_equality.matrix
		obj += call_f(phi)
	end

	for phi in dual_model.phi_equality.scalar
		obj += call_f(y, phi)
	end	

	obj += call_f(y, dual_model.phi_0)
	obj += call_f(y, dual_model.linear_term)

	return obj

end






















function call_f!(y::Vector{Float64}, phi::Phi_matrix; project::Bool=false)

	if project
		phi.M = projection(evaluate_matrix(y, phi)) # avoid using Matrix ?
	else
		phi.M = evaluate_matrix(y, phi) 
	end

	return -0.5*norm(phi.M)^2

end

function call_f_and_g!(y::Vector{Float64}, grad::Vector{Float64}, phi::Phi_matrix; project::Bool=false)

	if project
		phi.M = projection(evaluate_matrix(y, phi))
	else
		phi.M = evaluate_matrix(y, phi) 
	end

	grad[phi.address] = -[dot(Symmetric(phi.A_j[i]), phi.M) for i in 1:length(phi.A_j)] # avoid Symmetric ?
	
	return -0.5*norm(phi.M)^2 

end


function call_dual_obj!(y::Vector{Float64}, dual_model::DualModel)

	obj = call_f!(y, dual_model.phi_moment, project=true)

	for phi in dual_model.phi_inequality.matrix
		obj += call_f!(y, phi, project=true)
	end

	for phi in dual_model.phi_inequality.scalar
		obj += call_f(y, phi)
	end

	for phi in dual_model.phi_equality.matrix
		obj += call_f!(y, phi)
	end

	for phi in dual_model.phi_equality.scalar
		obj += call_f!(y, phi)
	end	

	obj += call_f(y, dual_model.phi_0)
	obj += call_f(y, dual_model.linear_term)

	return obj	 

end

function call_dual_obj_and_grad!(y::Vector{Float64}, grad::Vector{Float64}, dual_model::DualModel)

	obj, grad[:] = call_f_and_g!(y, grad, dual_model.phi_moment, project=true)

	for phi in dual_model.phi_inequality.matrix
		obj += call_f_and_g!(y, grad, phi, project=true)
	end

	for phi in dual_model.phi_inequality.scalar
		obj += call_f_and_g!(y, grad, phi)
	end

	for phi in dual_model.phi_equality.matrix
		obj += call_f_and_g!(y, grad, phi)
	end

	for phi in dual_model.phi_equality.scalar
		obj += call_f_and_g!(y, grad, phi)
	end	

	obj += call_f_and_g!(y, grad, dual_model.phi_0)
	obj += call_f_and_g!(y, grad, dual_model.linear_term)

	return obj

end
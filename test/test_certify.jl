using Test
using MomentHierarchy
MH = MomentHierarchy
using DynamicPolynomials

using LinearAlgebra
using SparseArrays

absapprox(x, y) = all(abs.(x - y) .<= 1e-7)

@testset "certification" begin

	n = 3
	@polyvar x[1:n]
	f = sum(x[i] for i in 1:n)
	g_pos = [x[i] + 0. for i in 1:n]
	g_sphere = sum(x[i]^2 for i in 1:n) - 1.
	pop = MH.POP(f, x, g_inequality=g_pos, g_equality=g_sphere)

	relaxation_order = 1
	val = 0.

	# phi moment matrix

	phi_moment, moment_labels = MH.set_phi_moment(pop, relaxation_order)

	@test length(phi_moment.A_j) == 10
	@test phi_moment.A_j[1] == dropzeros(sparse([1, 4], [1, 4], [1., 0.]))
	@test phi_moment.A_j[5] == dropzeros(sparse([2, 4], [3, 4], [1., 0.]))

	y = zeros(10)
	y[5] = 1.

	@test MH.evaluate_matrix(y, phi_moment) == Symmetric(dropzeros(sparse([2, 4], [3, 4], [1., 0.])))
	@test MH.projection(Symmetric(diagm([1., 2., -3.]))) == Symmetric(diagm([1., 2., 0.]))



	@test absapprox(MH.call_f!(y, phi_moment, project=false), -1.)
	@test absapprox(MH.call_f!(y, phi_moment, project=true), -0.5)

	grad = zeros(10)
	_ = MH.call_f_and_g!(y, grad, phi_moment, project=false)
	test_grad = zeros(10)
	test_grad[5] = -2.

	@test grad == test_grad

	grad = zeros(10)
	_ = MH.call_f_and_g!(y, grad, phi_moment, project=true)
	test_grad = zeros(10)
	test_grad[3] = -0.5
	test_grad[5] = -1.
	test_grad[6] = -0.5

	@test absapprox(grad, test_grad)

	# phi scalar

	@test MH.get_linear_form(pop.equality_constraints[1], moment_labels) == ([-1., 1., 1., 1.], [1, 3, 6, 10])

	phi_equality = MH.set_phi_scalar_equality(pop.equality_constraints[1], pop, moment_labels)

	y = zeros(10)
	y[10] = -2.

	@test MH.call_f(y, phi_equality) == -2.

	grad = zeros(10)
	_ = MH.call_f_and_g!(y, grad, phi_equality)
	test_grad = zeros(10)
	test_grad[1] = -2.
	test_grad[3] = 2.
	test_grad[6] = 2.
	test_grad[10] = 2.

	@test grad == test_grad

	phi_inequality = MH.set_phi_scalar_inequality(pop.inequality_constraints[1], pop, moment_labels)

	y = zeros(10)
	y[2] = -2.

	@test MH.call_f(y, phi_inequality) == 0.

	y[2] = 2.

	@test MH.call_f(y, phi_inequality) == -2.

	grad = zeros(10)
	_ = MH.call_f_and_g!(y, grad, phi_inequality)
	test_grad = zeros(10)
	test_grad[2] = -2.

	@test grad == test_grad

	# phi_0

	phi_0 = MH.set_phi_0(val)

	y[1] = 1.
	_ = MH.call_f_and_g!(y, grad, phi_0)

	@test MH.call_f(y, phi_0) == -0.5
	@test grad[1] == -1.

	y[1] = -1.
	grad[1] = 1.
	_ = MH.call_f_and_g!(y, grad, phi_0)
	
	@test MH.call_f(y, phi_0) == 0.
	@test grad[1] == 1.

	# linear term

	phi_linear_term = MH.set_linear_term(pop, moment_labels)

	y = zeros(10)
	y[7] = 2.
	grad = zeros(10)
	_ = MH.call_f_and_g!(y, grad, phi_linear_term)

	@test MH.call_f(y, phi_linear_term) == 2.
	@test grad[[2, 4, 7]] == ones(3) 

	dual_model = MH.build_dual_model(pop, relaxation_order, val)

	@test length(dual_model.phi_equality.scalar) == 1
	@test length(dual_model.phi_equality.matrix) == 0
	@test length(dual_model.phi_inequality.scalar) == 3
	@test length(dual_model.phi_inequality.matrix) == 0

	# phi localization matrix ### add test to check if equalities and inequalities are correct

	relaxation_order = 2

	phi_moment, moment_labels = MH.set_phi_moment(pop, relaxation_order)
	phi_inequality = MH.set_phi_matrix(pop.inequality_constraints[1], pop, 1, moment_labels)

	@test phi_inequality.A_j[2] == dropzeros(sparse([1, 4], [2, 4], [1., 0.]))
	@test phi_inequality.address[2] == 3

	p = MH.SparsePolynomial([[0x0001, 0x0001], [0x0001, 0x0001], [0x0000], [0x0000]], [1., 1., 1., 1.])
	phi_inequality = MH.set_phi_matrix(p, pop, 1, moment_labels)

	# ???
	@test phi_inequality.A_j[1] == dropzeros(sparse([1, 4], [1, 4], [2., 0.]))
	@test phi_inequality.address[1] == 1
	@test phi_inequality.A_j[2] == dropzeros(sparse([1, 4], [2, 4], [2., 0.]))
	@test phi_inequality.address[2] == 2


	@test sparse([5, 1, 1], [5, 1, 1], [1, 1, 1]) == sparse([5, 1], [5, 1], [1, 2])

	# add test on localization matrix for polynomial with redundants monomials (ex: 2*x[1] + 3*x[1])


end
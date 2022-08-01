# julia 1.6

using Test
using MomentHierarchy
MH = MomentHierarchy
using DynamicPolynomials
using JuMP

@testset "polynomials" begin
	
	@polyvar x[1:3]
	f = x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3]
	
	sparse_polynomial = MH.SparsePolynomial(f, x)

	@test sparse_polynomial.coefficients == [1.; 1.; 1.; 1.]
	@test sparse_polynomial.support == [[0x0001, 0x0001, 0x0001, 0x0001], 
										[0x0002, 0x0002, 0x0002, 0x0002], 
										[0x0003, 0x0003, 0x0003, 0x0003], 
										[0x0001, 0x0002, 0x0003]]
	@test MH.degree(sparse_polynomial) == 4

	pop = MH.POP(f, x, g_inequality=[- x[1]^2 - 0.5*x[2]^2, - x[2]^2 - x[3]^2])

	@test pop.n_variables == 3
	@test pop.objective.support == sparse_polynomial.support
	@test pop.inequality_constraints[1].coefficients == [-1. ; -0.5]
	@test pop.equality_constraints == nothing
	@test collect(MH.terms(sparse_polynomial)) == [([0x0001, 0x0001, 0x0001, 0x0001], 1.),
												   ([0x0002, 0x0002, 0x0002, 0x0002], 1.),
												   ([0x0003, 0x0003, 0x0003, 0x0003], 1.),
												   ([0x0001, 0x0002, 0x0003], 1.)]

end

@testset "moments" begin
	
	@polyvar x[1:3]
	f = x[1]^4 + 1
	sparse_f = MH.SparsePolynomial(f, x)
	pop = MH.POP(f, x)

	@test MH.n_moments(3, 1) == 4
	@test MH.n_moments(collect(1:3), 1) == 4
	@test MH.n_moments(pop, 2) == 10
	@test MH.localizing_matrix_order(sparse_f, 2) == 0

	@test collect(MH.moment_columns([0x0001], 2)) == [(1, [0x0000, 0x0000]),
													  (2, [0x0000, 0x0001]),
													  (3, [0x0001, 0x0001])]

	@test collect(MH.moment_rows([0x0001], 1, 2)) == [(1, [0x0000]),
 													  (2, [0x0001])]

 	@test collect(MH.coupling_moments([0x0001, 0x0003], 2)) == [[0x0000, 0x0000],
																[0x0000, 0x0001],
																[0x0000, 0x0003],
																[0x0001, 0x0001],
																[0x0001, 0x0003],
																[0x0003, 0x0003]]

	@test MH.monomial([0x0002, 0x0000, 0x0001]) == [0x0001, 0x0002]
	@test MH.monomial_product([0x0001], [0x0001, 0x0000]) == [0x0001, 0x0001]
	@test MH.monomial_product([0x0001], [0x0001, 0x0000], [0x0000]) == [0x0001, 0x0001]

end

@testset "dense model" begin
	
	@polyvar x[1:6]
	f = x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6]
	g_1 = 1 - x[1]^2 - x[2]^2 
	g_2 = 1 - x[3]^2 - x[4]^2  
	g_3 = 1 - x[5]^2 - x[6]^2 
	pop = MH.POP(f, x, g_inequality=g_1, g_equality=[g_2, g_3])
	
	relaxation_order = 1
	model = Model()

	MH.set_moment_variables!(model, pop, relaxation_order)
	@test length(model[:y]) == 28
	
	moment_labels = MH.set_moment_matrix!(model, pop, relaxation_order)
	@test length(keys(moment_labels)) == 28
	
	@test MH.linear_expression(model, pop.objective, moment_labels) == model[:y][5] + 
	model[:y][9] + model[:y][14] + model[:y][20] + model[:y][27]
	
	MH.set_objective!(model, pop, moment_labels)
	@test objective_function(model) == model[:y][5] + model[:y][9] + 
		model[:y][14] + model[:y][20] + model[:y][27]

	MH.set_probability_measure_constraint!(model, moment_labels)
	MH.set_polynomial_constraints!(model, pop, relaxation_order, moment_labels)

	@test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == 3
	@test num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) == 1
	@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 1

end

@testset "decomposition tools" begin
	
	@polyvar x[1:6]
	f = x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6]
	g = 1 - x[1]^2 - x[2]^2
	pop = MH.POP(f, x, g_equality=g)
	variable_sets = [[1, 2, 3], [3, 4, 5, 6]]

	relaxation_order = 1
	model = Model()

	MH.set_moment_variables!(model, variable_sets[1], 1, relaxation_order)
	@test length(model[:y_1]) == 10
	
	moment_labels = MH.set_moment_matrix!(model, variable_sets[1], 1, relaxation_order)
	@test length(keys(moment_labels)) == 10

	@test MH.objective_decomposition(pop, variable_sets)[1].support == [[0x0001, 0x0002],
																  		[0x0002, 0x0003]]
	@test MH.assign_constraint_to_set(pop.equality_constraints[1], variable_sets) == 1

end

@testset "dual decomposition" begin
	
	@testset "subproblems" begin



	end

	@testset "master" begin

	end


end

@testset "dummy decomposition" begin
	
	


end

@testset "nlp" begin
	
	
	
end
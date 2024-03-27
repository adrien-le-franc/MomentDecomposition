# julia 1.6

using Test
using MomentSOS
MSOS = MomentSOS
using DynamicPolynomials
using JuMP
using LinearAlgebra

@testset "basics" begin

	@testset "polynomials" begin
		
		@polyvar x[1:3]
		f = x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3]
		
		sparse_polynomial = MSOS.SparsePolynomial(f, x)

		@test sparse_polynomial.coefficients == [1., 1., 1., 1.]
		@test sparse_polynomial.support == [[0x0001, 0x0001, 0x0001, 0x0001], 
											[0x0002, 0x0002, 0x0002, 0x0002], 
											[0x0003, 0x0003, 0x0003, 0x0003], 
											[0x0001, 0x0002, 0x0003]]
		@test MSOS.degree(sparse_polynomial) == 4

		pop = MSOS.POP(f, x, g_inequality=[- x[1]^2 - 0.5*x[2]^2, - x[2]^2 - x[3]^2])

		@test pop.n_variables == 3
		@test pop.objective.support == sparse_polynomial.support
		@test pop.inequality_constraints[1].coefficients == [-1. ; -0.5]
		@test pop.equality_constraints == nothing
		@test collect(MSOS.terms(sparse_polynomial)) == [([0x0001, 0x0001, 0x0001, 0x0001], 1.),
													   ([0x0002, 0x0002, 0x0002, 0x0002], 1.),
													   ([0x0003, 0x0003, 0x0003, 0x0003], 1.),
													   ([0x0001, 0x0002, 0x0003], 1.)]

	end

	@testset "moments" begin
		
		@polyvar x[1:3]
		f = x[1]^4 + 1
		sparse_f = MSOS.SparsePolynomial(f, x)
		pop = MSOS.POP(f, x)

		@test MSOS.n_moments(3, 1) == 4
		@test MSOS.n_moments(collect(1:3), 1) == 4
		@test MSOS.n_moments(pop, 2) == 10
		@test MSOS.localizing_matrix_order(sparse_f, 2) == 0

		@test collect(MSOS.moment_columns([0x0001], 2)) == [(1, [0x0000, 0x0000]),
														  (2, [0x0000, 0x0001]),
														  (3, [0x0001, 0x0001])]

		@test collect(MSOS.moment_rows([0x0001], 1, 2)) == [(1, [0x0000]), (2, [0x0001])]
		@test MSOS.monomial([0x0002, 0x0000, 0x0001]) == [0x0001, 0x0002]
		@test MSOS.monomial_product([0x0001], [0x0001, 0x0000]) == [0x0001, 0x0001]
		@test MSOS.monomial_product([0x0001], [0x0001, 0x0000], [0x0000]) == [0x0001, 0x0001]

	end

end

@testset "models" begin
	
	@polyvar x[1:6]
	f = x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6]
	g_1 = 1 - x[1]^2 - x[2]^2 
	g_2 = 1 - x[3]^2 - x[4]^2  
	g_3 = 1 - x[5]^2 - x[6]^2 
	pop = MSOS.POP(f, x, g_inequality=g_1, g_equality=[g_2, g_3])
	
	relaxation_order = 1

	@testset "moment" begin
		
		model = Model()
		variable_sets = [collect(1:6)]
		
		moment_labels = MSOS.set_moment_matrices!(model, pop, relaxation_order, variable_sets)
		@test length(keys(moment_labels)) == 28
		
		@test MSOS.linear_expression(model, pop.objective, moment_labels) == model[:y][5] + 
		model[:y][9] + model[:y][14] + model[:y][20] + model[:y][27]
		
		MSOS.set_objective!(model, pop, moment_labels)
		@test objective_function(model) == model[:y][5] + model[:y][9] + 
			model[:y][14] + model[:y][20] + model[:y][27]

		MSOS.set_probability_measure_constraint!(model, moment_labels)
		MSOS.set_polynomial_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)

		@test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == 3
		@test num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) == 1
		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 1

		model = Model()
		moment_labels = MSOS.set_moment_matrices!(model, pop, 2, variable_sets)
		MSOS.set_polynomial_constraints!(model, pop, 2, moment_labels, variable_sets)

		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2	

	end

	@testset "sos" begin
		
		model = Model()
		@variable(model, t)
		variable_sets = [collect(1:6)]
		
		monomial_index = MSOS.set_X_0!(model, variable_sets, relaxation_order)
		@test length(keys(monomial_index)) == 28
		@test length(model[:X_0]) == 1

		MSOS.set_X_j!(model, pop, monomial_index, variable_sets, relaxation_order)
		@test length(model[:X_ineq]) == 1
		@test length(model[:X_eq]) == 2

		MSOS.set_SOS_constraints!(model, pop, monomial_index)
		@test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == 28
		@test num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) == 1
		@test num_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle) == 1

	end	

	@testset "nlp" begin
		
		model = MSOS.non_linear_model(pop)
		@test num_constraints(model; count_variable_in_set_constraints=false) == 3

	end

	@testset "certification" begin

		X = Symmetric([1 0 0; 0 1 0; 0 0 -2.0])
		@test MSOS.project_to_PSD_cone(X) == Symmetric([1 0 0; 0 1 0; 0 0 0.0])

		@polyvar x[1:2]
		f = x[1] + 0.
		g_1 = 1 - x[1]^2 - x[2]^2 
		pop = MSOS.POP(f, x, g_equality=g_1)
		model = MSOS.sos_relaxation_model(pop, 1)

		values = Dict{VariableRef, Float64}()

		values[model[:t]] = -1.0

		values[model[:X_0][1][1, 1]] = 1.0
		values[model[:X_0][1][1, 2]] = 0.5
		values[model[:X_0][1][1, 3]] = 0.
		values[model[:X_0][1][2, 2]] = 0.
		values[model[:X_0][1][2, 3]] = 1.0e-6
		values[model[:X_0][1][3, 3]] = 0.

		values[model[:X_eq][1]] = 0.

		monomial_index = Dict([0x0001] => 2)

		@test MSOS.evaluate_linear_system_error(model, pop, monomial_index, values) == -2.0e-6

	end

end

@testset "sparsity" begin

		@polyvar x[1:6]
		f = x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6]
		g_1 = 1 - x[1]^2 - x[2]^2 
		g_2 = 1 - x[3]^2 - x[4]^2  
		g_3 = 1 - x[5]^2 - x[6]^2 
		pop = MSOS.POP(f, x, g_inequality=g_1, g_equality=[g_2, g_3])

		model = Model()
		variable_sets = [[1, 2, 3], [3, 4, 5, 6]]
		relaxation_order = 1

		moment_labels = MSOS.set_moment_matrices!(model, pop, relaxation_order, variable_sets)
		
		@test length(keys(moment_labels)) == 22
		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2

		@test MSOS.assign_constraint_to_set(pop.inequality_constraints[1], variable_sets) == 1

	end
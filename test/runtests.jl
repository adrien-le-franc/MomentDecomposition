# julia 1.6

using Test
using MomentHierarchy
MH = MomentHierarchy
using DynamicPolynomials
using JuMP
using LinearAlgebra

@testset "basics" begin

	@testset "polynomials" begin
		
		@polyvar x[1:3]
		f = x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3]
		
		sparse_polynomial = MH.SparsePolynomial(f, x)

		@test sparse_polynomial.coefficients == [1., 1., 1., 1.]
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

end

@testset "models" begin
	
	@polyvar x[1:6]
	f = x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6]
	g_1 = 1 - x[1]^2 - x[2]^2 
	g_2 = 1 - x[3]^2 - x[4]^2  
	g_3 = 1 - x[5]^2 - x[6]^2 
	pop = MH.POP(f, x, g_inequality=g_1, g_equality=[g_2, g_3])
	
	relaxation_order = 1

	@testset "dense model" begin
		
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

		model = Model()
		MH.set_moment_variables!(model, pop, 2)
		moment_labels = MH.set_moment_matrix!(model, pop, 2)
		MH.set_polynomial_constraints!(model, pop, 2, moment_labels)

		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2	

	end

	variable_sets = [[1, 2, 3], [3, 4, 5, 6]]

	@testset "decomposition tools" begin
		
		@test MH.assign_constraint_to_set(pop.inequality_constraints[1], variable_sets) == 1
		@test MH.assign_constraint_to_sets(pop.inequality_constraints[1], variable_sets) == [1]
		@test collect(MH.pairs(variable_sets)) == [[1, 2]]
		@test MH.intersect(variable_sets, [1, 2]) == [3]

	end

	@testset "sparsity" begin

		model = Model()
		variable_sets = [[1, 2, 3], [3, 4, 5, 6]]

		MH.set_moment_variables!(model, pop, relaxation_order)
		moment_labels = MH.set_moment_matrices!(model, pop, relaxation_order, variable_sets)
		
		@test length(keys(moment_labels)) == 22
		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2

	end

	@testset "dual decomposition" begin
		
		@testset "subproblems" begin

			subproblems = [MH.SubProblem(k) for k in 1:2] 
			moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:2)

			for (k, set) in enumerate(variable_sets)

				MH.set_moment_variables!(subproblems[k].model, set, relaxation_order)
				moment_labels[k] = MH.set_moment_matrix!(subproblems[k].model, set, relaxation_order)
				MH.set_probability_measure_constraint!(subproblems[k].model, moment_labels[k]) 

			end

			@test length(subproblems[1].model[:y]) == 10
			@test length(moment_labels[1]) == 10

			MH.set_polynomial_objective!(subproblems, pop, moment_labels, variable_sets)
			p = AffExpr(0.)
			add_to_expression!(p, subproblems[1].model[:y][5] + subproblems[1].model[:y][9])
			@test subproblems[1].polynomial_objective == p

			MH.set_polynomial_constraints!(subproblems, pop, relaxation_order, # test order 2 ?
				moment_labels, variable_sets)

			@test num_constraints(subproblems[1].model, AffExpr, MOI.EqualTo{Float64}) == 1 
			@test num_constraints(subproblems[1].model, AffExpr, MOI.GreaterThan{Float64}) == 1
			@test num_constraints(subproblems[1].model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 1 

			info = MH.set_coupling_terms!(subproblems, [1, 2], 7, [3], 1, moment_labels, 2)

			@test subproblems[1].coupling_terms[1][1] == 1*subproblems[1].model[:y][7]
			@test subproblems[2].coupling_terms[1][1] == -1*subproblems[2].model[:y][2]
			@test length(subproblems[1].coupling_terms[1]) == 2
			@test subproblems[1].multiplier_ids == [7]
			@test info.moments == [[0x0003], [0x0003, 0x0003]]
			@test info.pair == (1, 2)

			subproblems = [MH.SubProblem(k) for k in 1:2] 
			moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:2)

			for (k, set) in enumerate(variable_sets)

				MH.set_moment_variables!(subproblems[k].model, set, 2)
				moment_labels[k] = MH.set_moment_matrix!(subproblems[k].model, set, 2)
				MH.set_probability_measure_constraint!(subproblems[k].model, moment_labels[k]) 

			end

			MH.set_polynomial_constraints!(subproblems, pop, 2,	moment_labels, variable_sets)

			@test num_constraints(subproblems[1].model, AffExpr, MOI.GreaterThan{Float64}) == 0
			@test num_constraints(subproblems[1].model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2

		end

		@testset "master" begin

			subproblems, info = MH.dual_decomposition_subproblems(pop, 1, variable_sets, max_coupling_order=1)
			@test info[1].moments == [[0x0003]]
			@test info[1].pair == (1, 2)

		end

		@testset "oracle" begin

			subproblems, info = MH.dual_decomposition_subproblems(pop, 1, variable_sets)

			multiplier = MH.Multiplier([[1.2, 3.4]], info)
			
			@test MH.extract(multiplier, subproblems[1]) == [[1.2, 3.4]]

			scale_factor = MH.update_dual_objective!(subproblems[1], [[2., -2.]])
			@test scale_factor == 2.
			@test objective_function(subproblems[1].model) == 0.5*subproblems[1].model[:y][5] + 
				0.5*subproblems[1].model[:y][9] + 1.0*subproblems[1].model[:y][7] - 1.0*subproblems[1].model[:y][10]

			oracle_data = [ [[1., 2.], [0.2]], [[3., 4.], [0.3]] ]

			@test MH.dual_objective_value(oracle_data) == 0.5
			@test MH.supergradient(subproblems, multiplier, oracle_data) == [[4., 6.]]
			@test MH.supergradient_norm(subproblems, multiplier, oracle_data) == norm([4., 6.])

		end

	end

	@testset "dummy decomposition" begin
		
		model = Model()
		moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:2)

		MH.set_moment_variables!(model, variable_sets[1], 1, relaxation_order)
		MH.set_moment_variables!(model, variable_sets[2], 2, relaxation_order)
		@test length(model[:y_1]) == 10
		
		moment_labels[1] = MH.set_moment_matrix!(model, variable_sets[1], 1, relaxation_order)
		moment_labels[2] = MH.set_moment_matrix!(model, variable_sets[2], 2, relaxation_order)
		@test length(keys(moment_labels[1])) == 10
		
		@test MH.linear_expression(model, pop.inequality_constraints[1], moment_labels, 1) == model[:y_1][1] - 
			model[:y_1][3] - model[:y_1][6] #
		
		MH.set_objective!(model, pop, moment_labels, variable_sets)

		@test objective_function(model) == model[:y_1][5] + model[:y_1][9] + 
			model[:y_2][5] + model[:y_2][9] + model[:y_2][14]

		MH.set_probability_measure_constraint!(model, moment_labels, 1)
		MH.set_probability_measure_constraint!(model, moment_labels, 2)

		MH.set_polynomial_constraints!(model, pop, relaxation_order, moment_labels, variable_sets)

		MH.set_coupling_constraints!(model, relaxation_order, moment_labels, variable_sets, 1)

		@test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == 2 + 2 + 1
		@test num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) == 1
		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2

		model = Model()
		moment_labels = Dict(k => Dict{Vector{UInt16}, Int64}() for k in 1:2)
		MH.set_moment_variables!(model, variable_sets[1], 1, 2)
		MH.set_moment_variables!(model, variable_sets[2], 2, 2)
		moment_labels[1] = MH.set_moment_matrix!(model, variable_sets[1], 1, 2)
		moment_labels[2] = MH.set_moment_matrix!(model, variable_sets[2], 2, 2)
		MH.set_polynomial_constraints!(model, pop, 2, moment_labels, variable_sets)

		@test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 3

	end

	@testset "nlp" begin
		
		model = MH.non_linear_model(pop)
		@test num_constraints(model; count_variable_in_set_constraints=false) == 3
		
	end

end
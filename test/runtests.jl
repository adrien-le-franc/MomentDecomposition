# julia 1.6

using Test
using MomentDecomposition
MD = MomentDecomposition
using DynamicPolynomials

@testset "polynomials" begin
	
	@polyvar x[1:3]
	f = x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3]
	
	sparse_polynomial = MD.SparsePolynomial(f, x)

	@test sparse_polynomial.coefficients == [1.; 1.; 1.; 1.]
	@test sparse_polynomial.support == [convert.(UInt16, [1; 1; 1; 1]), 
										convert.(UInt16, [2; 2; 2; 2]), 
										convert.(UInt16, [3; 3; 3; 3]), 
										convert.(UInt16, [1; 2; 3])] 

	pop = MD.POP(f, x, g_inequality=[- x[1]^2 - 0.5*x[2]^2, - x[2]^2 - x[3]^2])

	@test pop.n_variables == 3
	@test pop.objective.support == sparse_polynomial.support
	@test pop.inequality_constraints[1].coefficients == [-1. ; -0.5]
	@test pop.equality_constraints == nothing

end
# julia 1.6
#
# nonlinear model for polynomial optimization


function polynomial_expression(model::Model, polynomial::SparsePolynomial)

	support = [convert.(Int64, monomial) for monomial in polynomial.support]
	n_terms = length(support)

	expression = @NLexpression(model, 
		sum(polynomial.coefficients[i]*prod(model[:x][k] for k in support[i]) for i in 1:n_terms))

	return expression

end

function non_linear_model(pop::POP)

	model = Model()

	@variable(model, x[1:pop.n_variables])

	for g_inequality in pop.inequality_constraints
		polynomial = polynomial_expression(model, g_inequality)
		@NLconstraint(model, polynomial >= 0.)
	end

	for g_equality in pop.equality_constraints
		polynomial = polynomial_expression(model, g_equality)
		@NLconstraint(model, polynomial == 0.)
	end

	polynomial = polynomial_expression(model, pop.objective)
	@NLobjective(model, Min, polynomial)
	
	return model

end


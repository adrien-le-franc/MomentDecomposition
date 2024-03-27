# julia 1.6
#
# tools for correlative sparsity
#
# WARNING : when the relaxation order is 1 the csp graph 
# can be made even sparser than what build_csp_graph implements currently (to do)

function build_csp_graph(pop::POP)

	csp_graph = SimpleGraph(convert(UInt16, pop.n_variables))

	# objective

	for monomial in pop.objective.support

		variables = unique(monomial)

		for (i, j) in combinations(variables, 2)  

			add_edge!(csp_graph, i, j)		

		end

	end

	# constraints

	for polynomial in pop.inequality_constraints

		variables = unique(vcat(polynomial.support...))

		for (i, j) in combinations(variables, 2)  

			add_edge!(csp_graph, i, j)		

		end

	end

	for polynomial in pop.equality_constraints

		variables = unique(vcat(polynomial.support...))

		for (i, j) in combinations(variables, 2)  

			add_edge!(csp_graph, i, j)		

		end

	end

	return csp_graph

end

function compute_cliques(pop::POP)
	
	csp_graph = build_csp_graph(pop)
	cliques, num_cliques, size_cliques = chordal_cliques!(csp_graph, method="MF", minimize=true)
	return cliques

end
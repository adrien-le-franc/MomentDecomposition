# julia 1.6
#
# tools for term sparsity


function compute_tsp_blocks_X_0(sparsity_pattern::SparsityPattern, k::Int64, relaxation_order::Int64)

	basis_size = n_moments(sparsity_pattern.variable_sets[k], relaxation_order)
	tsp_graph = SimpleGraph(basis_size)

	for (j, alpha) in moment_columns(sparsity_pattern.variable_sets[k], relaxation_order)
		for (i, beta) in moment_rows(sparsity_pattern.variable_sets[k], relaxation_order, j)

			monomial = monomial_product(alpha, beta)

			if monomial in sparsity_pattern.monomial_sets[k] || all(iseven.(degrees(monomial)))
				if i != j
					add_edge!(tsp_graph, i, j)
				end
			end 

		end
	end

	return connected_components(tsp_graph)

end

function scalar(blocks) 
	return enumerate(@view(blocks[length.(blocks) .== 1]))	
end

function matrix(blocks)
	return enumerate(@view(blocks[length.(blocks) .> 1]))
end
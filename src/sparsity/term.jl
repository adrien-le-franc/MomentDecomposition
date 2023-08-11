# julia 1.6
#
# tools for term sparsity


function compute_tsp_blocks_X_0(variable_set, monomial_set, relaxation_order::Int64)

	basis_size = n_moments(variable_set, relaxation_order)
	tsp_graph = SimpleGraph(basis_size)

	for (j, alpha) in moment_columns(variable_set, relaxation_order)
		for (i, beta) in moment_rows(variable_set, relaxation_order, j)

			monomial = monomial_product(alpha, beta)

			if monomial in monomial_set || all(iseven.(degrees(monomial)))
				if i != j
					add_edge!(tsp_graph, i, j)
				end
			end 

		end
	end

	return connected_components(tsp_graph)

end

function diagonal_block(blocks)
	return Vector{Int64}(@view(blocks[length.(blocks) .== 1]))
end

function matrix_blocks(blocks)
	return enumerate(@view(blocks[length.(blocks) .> 1]))
end
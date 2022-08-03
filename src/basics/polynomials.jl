# julia 1.6
#
# tools for large scale polynomial optimization

struct SparsePolynomial
	support::Vector{Vector{UInt16}}
	coefficients::Vector{Float64}
end

function SparsePolynomial(f::P, variables::Vector{V}) where {P <: Polynomial, V <: PolyVar}

	n = length(variables)
    
    terms = DynamicPolynomials.monomials(f)
    n_temrs = length(terms)
    coefficients = DynamicPolynomials.coefficients(f)
    support = [UInt16[] for i=1:n_temrs]

    for i in 1:n_temrs
        indices = terms[i].z .> 0
        term_variables = terms[i].vars[indices]
        exponents = terms[i].z[indices]
        for j in 1:length(term_variables)
            variable_indice = find_variable_indice(variables, n, term_variables[j])
            append!(support[i], variable_indice*ones(UInt16, exponents[j]))
        end
    end

    return SparsePolynomial(support, coefficients) 

end

degree(f::SparsePolynomial) = maximum(length(monomial) for monomial in f.support)
terms(f::SparsePolynomial) = zip(f.support, f.coefficients)

struct POP
	objective::SparsePolynomial
	n_variables::Int64
	inequality_constraints::Union{Nothing, Vector{SparsePolynomial}}
	equality_constraints::Union{Nothing, Vector{SparsePolynomial}}
end 

function POP(f::P1, x::Vector{V};
	g_inequality::Union{Nothing, P2, Vector{P2}}=nothing,
	g_equality::Union{Nothing, P3, Vector{P3}}=nothing) where {P1 <: Polynomial, V <: PolyVar,
		P2 <: Polynomial, P3 <: Polynomial}
	
	objective = SparsePolynomial(f, x)
	n_variables = length(x)

	if typeof(g_inequality) <: Vector 
		inequality_constraints = [SparsePolynomial(g, x) for g in g_inequality]
	elseif typeof(g_inequality) <: Polynomial
		inequality_constraints = [SparsePolynomial(g_inequality, x)]
	else
		inequality_constraints = nothing
	end

	if typeof(g_equality) <: Vector 
		equality_constraints = [SparsePolynomial(g, x) for g in g_equality]
	elseif typeof(g_equality) <: Polynomial
		equality_constraints = [SparsePolynomial(g_equality, x)]
	else
		equality_constraints = nothing
	end

	return POP(objective, n_variables, inequality_constraints, equality_constraints)

end

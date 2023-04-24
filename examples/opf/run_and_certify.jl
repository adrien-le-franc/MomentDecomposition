# julia 1.6

ENV["MOSEKLM_LICENSE_FILE"] = "/home/OPF/mosek/mosek.lic"
#DATA_FOLDER = "/home/OPF/data/pglib-opf/"
DATA_FOLDER = "/home/OPF/data/case57_difficult/case57/"
RESULT_FOLDER = "/home/OPF/MomentHierarchy.jl/x/opf/results"
max_set_size = 12
relaxation_order = 2
case57_rte = true
#case = "pglib_opf_case5_pjm.m"
case = "case57_391.m"


using MomentSOS
MSOS = MomentSOS

using MosekTools
using Ipopt
using JuMP

using FileIO
using Dates

include("certify.jl")
include("scaling.jl")
include("parse.jl")
include("parse_linear.jl")

# AC-OPF to POP

data = parse_file(joinpath(DATA_FOLDER, case))

if case57_rte
	pop, x0, minimal_sets, bounds = parse_opf_linear_costs_to_pop(data)
else
	pop, x0, minimal_sets, bounds = parse_opf_to_pop(data, AngleCons=true, LineLimit=true, nlp=true, n_max=max_set_size)
end

objective_scaling, objective = scale_polynomial(pop.objective, bounds, true)
inequalities = [scale_polynomial(g, bounds, true)[2] for g in pop.inequality_constraints]
equalities = [scale_polynomial(g, bounds, true)[2] for g in pop.equality_constraints]
pop = MSOS.POP(objective, pop.n_variables, inequalities, equalities)

# NLP

model = MSOS.non_linear_model(pop)
set_optimizer(model, Ipopt.Optimizer)
set_start_value.(model[:x], x0)
optimize!(model)

nlp = Dict("primal_status" => primal_status(model),
		   "dual_status" => dual_status(model),
		   "termination_status" => termination_status(model),
		   "upper_bound" => objective_value(model)*objective_scaling)

# Relaxation 

t_build = @elapsed model, monomial_index = MSOS.sos_relaxation_model(pop, relaxation_order, minimal_sets, true)

set_optimizer(model, Mosek.Optimizer)
optimize!(model)
t_solve = solve_time(model)

t_error = @elapsed err = compute_error(model, pop, monomial_index)
lower_epsilon = sum([e <= 0. ? e : 0. for e in err])
certified_bound = (objective_value(model) + lower_epsilon)*objective_scaling


relaxation = Dict("t_build" => t_build, "t_solve" => t_solve,
				  "primal_status" => primal_status(model),
				  "dual_status" => dual_status(model),
				  "termination_status" => termination_status(model),
				  "lower_bound" => objective_value(model)*objective_scaling,
				  "certified_bound" => certified_bound)

# Results

println("NLP")
display(nlp)
println("Relaxation")
display(relaxation)
# dense hierarchy resolution for 
#
# min f(x) st x in K 
#  x in |R^2
#
# where f = x[1]^3 - x[1]^2 + 2*x[1]*x[2] - x[2]^2 + x[2]^3
#
# and K is defined by:
#
# g_1 = x[1] >= 0
# g_2 = x[2] >= 0
# g_3 = x[1]^3 - 4*x[1]*x[2] + x[2]^2 - 1. == 0 


using DynamicPolynomials
using JuMP
using MomentSOS
MSOS = MomentSOS

ENV["MOSEKLM_LICENSE_FILE"] = "/home/OPF/mosek/mosek.lic"
using MosekTools

@polyvar x[1:5]

f = x[1]^3 - x[1]^2 + 2*x[1]*x[2] - x[2]^2 + x[2]^3 + x[3]^2 + x[4]^2 + 7*x[5]^2

g_1 = x[1] + 0.
g_2 = x[2] + 0.
g_3 = 4*x[3] - x[4]^2 + x[5] + 1.

g_4 = x[1]^3 - 4*x[1]*x[2] + x[2]^2 - 1.
g_5 = x[1]*x[2] - 2*x[4] + x[2]*x[5]^2

pop = MSOS.POP(f, x, g_inequality=[g_1, g_2, g_3], g_equality=[g_4, g_5])

model = MSOS.moment_relaxation_model(pop, 2)

set_optimizer(model, Mosek.Optimizer)

optimize!(model)
println(objective_value(model))

# v_pop = 0 
# -----------
# v_3 = 4.55e-8
# v_2 = -0.007955
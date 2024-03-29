# MomentSOS.jl

This package implements JuMP models for the moment-SOS hierarchy in polynomial optimization. It is designed for research purposes, with a specific focus on the following properties: 

* **flexibility** (you can make changes and access the underlying JuMP model easily),
*  **tests** (as much as possible, the code is tested),
* **fidelity to the mathematical model** (if you are familiar with the moment-SOS hierarchy, you should be able to understand this code).

Other packages with more advanced functionalities include [TSSOS](https://github.com/wangjie212/TSSOS), [SumsOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl), [MomentTools.jl](https://github.com/AlgebraicGeometricModeling/MomentTools.jl) and [GloptiPoly](https://homepages.laas.fr/henrion/software/gloptipoly/).

## An introductive example

Let's consider the following Polynomial Optimization Problem (POP)

$$\min_{x \in \mathbb{R}^3} \quad x_1x_2 + x_2x_3 \quad \text{ s.t. } \quad 1 - x_3^2 \geq 0, \ x_1^2 + x_2^2 - 1 = 0. $$

We can encode the data of the problem with [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) 

```julia
using DynamicPolynomials

@polyvar x[1:3]
f = x[1]*x[2] + x[2]*x[3]
g_1 = 1 - x[3]^2 
g_2 = x[1]^2 + x[2]^2 - 1.
```

To solve the problem with the moment-SOS hierarchy, we load MomentSOS, JuMP and any [supported](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) SDP solver (in this example we use Clarabel)

```julia
using MomentSOS
MSOS = MomentSOS
using JuMP
using Clarabel
```

Then, we define the POP and the first-order SOS relaxation

```julia
pop = MSOS.POP(f, x, g_inequality=g_1, g_equality=g_2)
model = MSOS.sos_relaxation_model(pop, 1)
```

Alternatively, we may choose to work with the moment relaxation

```julia
model = MSOS.moment_relaxation_model(pop, 1)
```

In both cases, `model` is a JuMP model: we can set a SDP solver and solve the relaxation as follows

```julia
set_optimizer(model, Clarabel.Optimizer)
optimize!(model)
```

## Further examples

Further insights on more advanced functionalities are provided in the [examples](https://github.com/adrien-le-franc/MomentSOS.jl/tree/main/examples) folder. In particular, the code used for 
the paper [Minimal Sparsity for Second-Order Moment-SOSRelaxations of the AC-OPF Problem](https://hal.science/hal-04110742v2/document) is available in the [opf](https://github.com/adrien-le-franc/MomentSOS.jl/tree/main/examples/opf) subfolder.

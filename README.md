# Virtual Element Solution with Gridap.jl (WIP)

### **Note**: 
- **This code works only for the lowest order virtual element. Currently the code has the Ritz projector built in. This is good, since for the lowest order VEM, the Ritz projector can be used for both stiffness and mass terms.**
- **This code heavily depends on the linear Lagrange `FESpace` machinery in Gridap.jl. The `FESpace` machinery is used for the geometrical part and for assembly, as the degrees of freedom of the lowest order VEM is identical to linear Lagrange `FESpace`.**

### Example 1:

An example problem using the script `vem_example_gridap.jl`

```julia
include("VEM-Julia.jl")

using Gridap: ∇

domain = (0,1,0,1)
u(x) = sin(π*x[1])*sin(π*x[2])

K(x) = x[1]^2 + x[2]^2
σ(x) = K(x)⋅(∇(u))(x)
f(x) = -(∇⋅σ)(x) + u(x)

partition = (10,10)
model = CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 4)

Vₕ = P1ConformingVESpace(model, K; dirichlet_tags="boundary"); # Argument K to define the scaling of the stability term
Vₕ⁰ = TrialVESpace(Vₕ, 0);

a(u,v)=∫(K*∇(u)⋅∇(v) + u*v)Qₕ
l(v) = ∫(f*v)Qₕ

op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)

uh = solve(op) ## Get the FEFunction

e = u - uh
err = √(sum(∫(e*e)Qₕ))
```

### Example 2:

Check rate of convergence using the script `ooc_vem_example.jl`:
  ```
 julia> err
5×1 Matrix{Float64}:
 0.040594555904984765
 0.010306757323680199
 0.002586435110433769
 0.0006472156606070523
 0.00016184181015375544

julia> log.(err[1:end-1] ./ err[2:end])./log.(2)
4-element Vector{Float64}:
 1.9776957536320143
 1.994553606200162
 1.9986465750207745
 1.9996621557003598
  ```
  
[Check out my blog for details](https://balaje.github.io/2021/09/29/Virtual-Elements.html)
  

## Virtual Element Solution with Gridap.jl (WIP)

**Note**: 
- **This code works only for the lowest order virtual element. Currently the code has the Ritz projector built in. This is good, since for the lowest order VEM, the Ritz projector can be used for both stiffness and mass terms.**
- **This code heavily depends on the linear Lagrange `FESpace` machinery in Gridap.jl**

An example problem:

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

## Description:

### FESpace and Element-wise contributions
- To build a lowest order VE space for a discrete model, we use the new data structure
  ```julia
  P1ConformingVESpace(::DiscreteModel, ::StabilityCoeffs; kwargs...):
    model::DiscreteModel
    Π∇
    stability_term
    linear_fespace::FESpace # FESpace(model, ReferenceFE(lagrangian, Float64, 1); kwargs...)
    stab_coeff::Function
  
  ```
  For the lowest order order case, the degrees of freedom is identical to the lowest order Lagrange finite element space. So for geometrical purposes, we can associate a linear `FESpace` to the new conforming `VESpace`. 
  In addition to that we have two new data structures `Π∇, stability_term` for the purpose of (cellwise) VEM projectors. You can find the details of the construction [here](http://arturo.imati.cnr.it/brezzi/papers/hitchhikers-preprint.pdf).
  
- Once we construct the `VESpace`, we can construct the stiffness matrices (local) and the load vector (local) using these functions cell-wise:
  ```julia
  function _generate_mat_contribs(a::Function, V::VESpace, acd, ind)
    area, centroid, diameter = acd[end]
    stab_coeff = get_stab_coeff(V)
    :
    # Some magic here ... #
    :
    coeff_matrix = stab_coeff(centroid)
    projector'*H*projector + (coeff_matrix)*stability_term # Final term
  end
  ```
  Here the `coeff_matrix` is a diagonal matrix of the stability function evaluated at the element centroid. This defines the scaling of the stability term. In the example, for the weak formulation
  ```julia
  ∫(K*∇(u)⋅∇(v) + u*v)
  ```
  the diffusion function `K` is used to define the scaling of the term, since the weak form scales `~ K`. The load vector similarly 
  ```julia
  function _generate_vec_contribs(l::Function, V::VESpace, acd, ind)
    a, c, d = acd[ind]
    :
    # More magic ... #
    :
    projector'*Fⱼ
  end
  ```
  where we don't have the stability term. The `magic` here is that `Gridap` can be used to compute the action of the bilinear form `a(u,v)` on the scaled monomials. This is because the quadrature rules on the elements are well defined and direct numerical integration is possible (Do take a look at how the implementation works in the code [here](https://github.com/Balaje/VEM-Julia-Gridap/blob/740985a078551cf92d9815c3e48c319d54e8d150/src/VESpace/AffineVEOperator.jl#L48) and [here](https://github.com/Balaje/VEM-Julia-Gridap/blob/740985a078551cf92d9815c3e48c319d54e8d150/src/VESpace/AffineVEOperator.jl#L83)). However, note that this fails for arbitrary order polygons. The idea is then to subdivide the polygons into triangles and dispatch `∫` appropriately (Future work).
  
- To obtain the cell-wise contributions, we map the above functions to each element in the mesh: We do this using the `lazy_map` functionality.
  ```julia
    matcontribs = lazy_map(i -> _generate_mat_contribs(a, trial, acd, i), 1:num_cells(mesh))
    loadcontribs = lazy_map(i -> _generate_vec_contribs(l, trial, acd, i), 1:num_cells(mesh))
  ```
  
## Assembly and solution
We then use the low-level assembly routine: `SparseMatrixAssembler` available in Gridap, but with the extra dispatch
  ```julia
  function SparseMatrixAssembler(trial::VESpace, test::VESpace)
    SparseMatrixAssembler(trial.linear_fespace, test.linear_fespace)
  end
  
  # Assemble the matrices
  σₖ = get_cell_dof_ids(trial.linear_fespace)
  matdata = ([matcontribs],[σₖ],[σₖ])
  vecdata = ([loadcontribs], [σₖ])
  assem = SparseMatrixAssembler(trial, test)
  matrix = assemble_matrix(assem, matdata)
  vector = assemble_vector(assem, vecdata)
  ```
  
We then solve the problem and interpolate the raw vector (VEM solution) to the underlying `FESpace`, since we cannot express the virtual element basis explicitly. We perform a very basic convergence analysis to see if the numbers obey the error estimate
  ```julia
  || u - uₕ ||₀ ≤ Ch²
  ```
  We run the rate of convergence script `ooc_vem_example.jl` to check:
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
  
  

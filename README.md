## Virtual Element Solution with Gridap.jl (WIP)

**Note**: 
- **This code works only for the lowest order virtual element. Currently the code has the Ritz projector built in. This is good, since for the lowest order VEM, the Ritz projector can be used for both stiffness and mass terms.**
- **This code heavily depends on the linear Lagrange `FESpace` machinery in Gridap.jl**

An example problem:

```julia
include("VEM-Julia.jl") # Load the functions

domain = (0,1,0,1)
f(x) = (2π^2+1)*sin(π*x[1])*sin(π*x[2])
u(x) = sin(π*x[1])*sin(π*x[2])

partition = (50,50)
model = simplexify(CartesianDiscreteModel(domain, partition))
Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 4)

Vₕ = P1ConformingVESpace(model, [0,1]; dirichlet_tags="boundary");
Vₕ⁰ = TrialVESpace(Vₕ, 0);

# Simple weak form. Varies in VEM case (with consistency and stability terms ...), but is simple for the lowest order. 
a(u,v)=∫(∇(u)⋅∇(v) + u*v)Qₕ
l(v) = ∫(f*v)Qₕ

# Solve
op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)
uh = solve(op) ## Get the FEFunction
```

Description:

- To build a lowest order VE space for a discrete model, we use the new data structure
  ```julia
  P1ConformingVESpace(::DiscreteModel, ::StabilityCoeffs; kwargs...):
    model::DiscreteModel
    Π∇
    stability_term
    linear_fespace::FESpace # FESpace(model, ReferenceFE(lagrangian, Float64, 1); kwargs...)
    stab_coeff
  
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
    projector'*H*projector + sum(area.^stab_coeff)*stability_term # Final term
  end
  ```
  Here the stability term is defined as sums of powers of `area(element) = |K|`. You have `1` in every mass terms and `0` in every stiffness terms (involves gradients). So for this type of weak form
  ```julia
  ∫(∇(u)⋅∇(v) + u*v)
  ```
  you have the sum of `|K|^0` and `|K|^1` for the stiffness and mass terms respectively, and hence `Vₕ.stab_coeff = [0,1]`. The load vector similarly 
  ```julia
  function _generate_vec_contribs(l::Function, V::VESpace, acd, ind)
    a, c, d = acd[ind]
    :
    # More magic ... #
    :
    projector'*Fⱼ
  end
  ```
This is possible because Gridap allows us to do something like 

```julia
mᵢ(x) = x[1]; cfᵢ = CellField(mᵢ, mesh)
mⱼ(x) = x[2]; cfⱼ = CellField(mⱼ, mesh)

H[i,j] = ∫(∇(cfᵢ)⋅∇(cfⱼ) + cfᵢ*cfⱼ)Qₕ
```

This can be used to compute the action of the bilinear form `a(u,v)` on the scaled monomials. This is because the quadrature rules on the elements are well defined and direct numerical integration is possible. However, note that this fails for arbitrary order polygons. The idea is then to subdivide the polygons into triangles and dispatch `∫` appropriately.
  
- To obtain the cell-wise contributions, we map the above functions to each element in the mesh: We do this using the `lazy_map` functionality.
  ```julia
    matcontribs = lazy_map(i -> _generate_mat_contribs(a, trial, acd, i), 1:num_cells(mesh))
    loadcontribs = lazy_map(i -> _generate_vec_contribs(l, trial, acd, i), 1:num_cells(mesh))
  ```
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
  
- We then solve the problem and interpolate the raw vector (VEM solution) to the underlying `FESpace`, since we cannot express the virtual element basis explicitly. We perform a very basic convergence analysis to see if the numbers obey the error estimate
  ```julia
  || u - uₕ || ≤ Ch²
  ```
  We run the rate of convergence script `ooc_vem_example.jl` to check:
  ```
  julia> err
  5×1 Matrix{Float64}:
   0.05044987832391834
   0.013124575219057735
   0.003315286192444438
   0.000830998335419876
   0.0002078864156536598

  julia> log.(err[1:end-1] ./ err[2:end])./log.(2)
  4-element Vector{Float64}:
   1.9425800595509743
   1.985065408628436
   1.9962159242854818
   1.9990500988135191
  ```
  
  

include("VEM-Julia.jl")

using Gridap: ∇

domain = (0,1,0,1)
u(x) = sin(π*x[1])*sin(π*x[2])

K(x) = 1
σ(x) = K(x)⋅(∇(u))(x)
f(x) = -(∇⋅σ)(x) + u(x)


partition = (40,40)
model = CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 4)

# Pass K into the FESpace for the stability term
# - This could be any function that scales appropriately as the bilinear form.
# - The closest function is actually K

Vₕ = P1ConformingVESpace(model, K; dirichlet_tags="boundary");
Vₕ⁰ = TrialVESpace(Vₕ, 0);

a(u,v)=∫(K*∇(u)⋅∇(v) + u*v)Qₕ
l(v) = ∫(f*v)Qₕ

op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)

uh = solve(op) ## Get the FEFunction

e = u - uh
err = √(sum(∫(e*e)Qₕ))

writevtk(get_triangulation(uh), "ve_solution", cellfields=["gh"=>uh])
writevtk(get_triangulation(uh), "exact_solution", cellfields=["gh"=>CellField(u, Ω)])

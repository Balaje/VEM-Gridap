include("VEM-Julia.jl")

using Gridap: ∇

domain = (0,1,0,1)
u(x) = sin(π*x[1])*sin(π*x[2])

K(x) = x[1]^2 + x[2]^2
σ(x) = K(x)⋅(∇(u))(x)
f(x) = -(∇⋅σ)(x) + u(x)


partition = (20,20)
model = simplexify(CartesianDiscreteModel(domain, partition))
Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 4)

Vₕ = P1ConformingVESpace(model, [0,1]; dirichlet_tags="boundary");
Vₕ⁰ = TrialVESpace(Vₕ, 0);

a(u,v)=∫(K*∇(u)⋅∇(v) + u*v)Qₕ
l(v) = ∫(f*v)Qₕ

op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)

uh = solve(op) ## Get the FEFunction

e = u - uh
err = √(sum(∫(e*e)Qₕ))

writevtk(get_triangulation(uh), "ve_solution", cellfields=["gh"=>uh])
writevtk(get_triangulation(uh), "exact_solution", cellfields=["gh"=>CellField(u, Ω)])

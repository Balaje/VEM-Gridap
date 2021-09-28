include("VEM-Julia.jl")

domain = (0,1,0,1)
f(x) = (2π^2+1)*sin(π*x[1])*sin(π*x[2])
u(x) = sin(π*x[1])*sin(π*x[2])

partition = (50,50)
model = CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 4)

Vₕ = P1ConformingVESpace(model, [0,1]; dirichlet_tags="boundary");
Vₕ⁰ = TrialVESpace(Vₕ, 0);

a(u,v)=∫(∇(u)⋅∇(v) + u*v)Qₕ
l(v) = ∫(f*v)Qₕ

op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)

uh = solve(op) ## Get the FEFunction

writevtk(get_triangulation(uh), "ve_solution", cellfields=["gh"=>uh])
writevtk(get_triangulation(uh), "exact_solution", cellfields=["gh"=>CellField(u, Ω)])

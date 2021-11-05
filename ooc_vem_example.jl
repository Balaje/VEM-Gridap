include("VEM-Julia.jl")

using Plots
using LaTeXStrings

domain = (0,1,0,1)
K = 2 # Diffusion constant
f(x) = (2π^2*K+1)*sin(π*x[1])*sin(π*x[2])
u(x) = sin(π*x[1])*sin(π*x[2])

partition = [(5,5), (10,10), (20,20), (40,40), (80,80)]
max_ref = 5
err = zeros(max_ref, 1)

for nref = 1:max_ref
  model = simplexify(CartesianDiscreteModel(domain, partition[nref]))
  Ω = Triangulation(model)
  Qₕ = CellQuadrature(Ω, 4)

  # Pass the scaling function x->K (diffusion constant) to modify the stability term
  Vₕ = P1ConformingVESpace(model, x->K; dirichlet_tags="boundary");
  Vₕ⁰ = TrialVESpace(Vₕ, 0);

  a(u,v)=∫(K*∇(u)⋅∇(v) + u*v)Qₕ
  l(v) = ∫(f*v)Qₕ

  op = AffineVEOperator(a, l, Vₕ⁰, Vₕ)

  uh = solve(op) ## Get the FEFunction

  e = u - uh

  err[nref] = √(sum(∫(e*e)Qₕ))

  print("Done level ", nref, "\n")

  # Uncomment to visulaize
  # writevtk(get_triangulation(uh), "ve_solution", cellfields=["gh"=>uh])
  # writevtk(get_triangulation(uh), "exact_solution", cellfields=["gh"=>CellField(u, Ω)])
end

N = [5,10,20,40,80];
h = 1 ./N;

plt = plot(log.(h), log.(err), marker=5, label="\$|| u - u_h ||_{0}\$");
plot!(log.(h), log.(h.^2))

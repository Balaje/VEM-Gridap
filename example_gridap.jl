using Gridap

domain = (0,1,0,1)
partition = (20,20)
reffe = ReferenceFE(lagrangian, Float64, 1)
model = CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)
dΩ = Measure(Ω,2)

V = FESpace(model, reffe)
a(u,v) = ∫(u*v)dΩ

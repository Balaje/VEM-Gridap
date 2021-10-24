"""
Declare a struct containing the VEM space. Contains
- projectors Π∇, Π0
- stability term
The projectors and stability term depend only on the geometry of the element.
"""
abstract type VESpace end

struct P1ConformingVESpace <: VESpace
  model::DiscreteModel
  Π∇
  stability_term
  linear_fespace::FESpace
  stab_coeff::Function
end

struct P2ConformingVESpace <: VESpace
  model::DiscreteModel
  Π∇
  Π0
  stability_term
  linear_fespace::FESpace
  stab_coeff::Function
end

function P1ConformingVESpace(model::DiscreteModel, stab_coeff; kwargs...)
  mesh = Triangulation(model)
  cell_verts = get_cell_coordinates(mesh)
  Π∇ = lazy_map(_generate_ritz_matrices, cell_verts, geo(mesh))
  stability_term = lazy_map(_generate_stabilising_term, Π∇, geo(mesh))
  linear_fespace = FESpace(model, ReferenceFE(lagrangian, Float64, 1); kwargs...)
  P1ConformingVESpace(model, Π∇, stability_term, linear_fespace, stab_coeff)
end

function P2ConformingVESpace(model::DiscreteModel, stab_coeff; kwargs...)
  mesh = Triangulation(model)
  cell_verts = get_cell_coordinates(mesh)
  Π∇ = lazy_map(_generate_ritz_matrices, cell_verts, geo(mesh))
  Π0 = lazy_map(_generate_l2_matrices, cell_verts, geo(mesh))
  stability_term = lazy_map(_generate_stabilising_term, Π∇, geo(mesh))

end

struct TrialVESpace <: VESpace
  space::P1ConformingVESpace
  linear_fespace::FESpace
end

function TrialVESpace(f::VESpace, object)
  fespace = f.linear_fespace
  TrialVESpace(f, TrialFESpace(fespace, object))
end

FESpaces.get_triangulation(V::P1ConformingVESpace) = get_triangulation(V.linear_fespace)
get_Π∇(V::P1ConformingVESpace) = V.Π∇
get_stability(V::P1ConformingVESpace) = V.stability_term
get_stab_coeff(V::P1ConformingVESpace) = V.stab_coeff

FESpaces.get_triangulation(V::TrialVESpace) = get_triangulation(V.space.linear_fespace)
get_Π∇(V::TrialVESpace) = V.space.Π∇
get_stability(V::TrialVESpace) = V.space.stability_term
get_stab_coeff(V::TrialVESpace) = V.space.stab_coeff

"""
Routine for computing G, D, B, H for one element
"""
function _generate_ritz_matrices(verts, acd)
  if(length(verts) == 4)
    verts = view(verts,[1,2,4,3]) # To make it a closed quad
  end
  nsides = length(verts)
  npolys = 3
  a,c,d = acd
  modwrap = (x,a) -> mod(x-1, a) + 1
  D = zeros(nsides, npolys)
  D[:,1].=1
  B=zeros(npolys, nsides)
  B[1,:] .= 1/nsides
  linear_polynomials=[[0,0],[1,0],[0,1]]
  vs = map(nc -> SVector(Tuple(nc)), verts)
  vs = reinterpret(reshape, Float64, vs)

  for v=1:nsides
    vert = vs[:,v]
    prev = vs[:,modwrap(v-1, nsides)]
    next = vs[:,modwrap(v+1, nsides)]
    normal = Point(next[2]-prev[2], prev[1]-next[1])
    for poly_id = 2:npolys
      poly_degree = linear_polynomials[poly_id]
      monomial_grad = poly_degree/d
      D[v, poly_id] = ((Point(vert) .- c) ⋅ poly_degree)/d
      B[poly_id, v] = 0.5*(monomial_grad ⋅ normal)
    end
  end
  B, D
end


function _generate_stabilising_term(Π∇, acd)
  B, D = Π∇
  projector = (B*D)\B
  nsides = size(B,2)
  (I(nsides) - D*projector)'*(I(nsides) - D*projector)
end

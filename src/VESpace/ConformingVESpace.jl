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
  stab_coeff
end


function P1ConformingVESpace(model::DiscreteModel, stab_coeff; kwargs...)
  mesh = Triangulation(model)
  verts = mesh.node_coords
  m = Broadcasting(Reindex(verts))
  cell_verts = lazy_map(m, mesh.cell_node_ids)

  Π∇ = lazy_map(_generate_ritz_matrices, cell_verts, geo(mesh))
  stability_term = lazy_map(_generate_stabilising_term, Π∇, geo(mesh))

  linear_fespace = FESpace(model, ReferenceFE(lagrangian, Float64, 1); kwargs...)
  P1ConformingVESpace(model, Π∇, stability_term, linear_fespace, stab_coeff)
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
  verts = view(verts,[1,2,4,3]) # To make it a closed polygon.
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

function _generate_l2_matrices(verts, acd)
  npolys = 3
  H = zeros(npolys, npolys)
  vs = map(nc -> SVector(Tuple(nc)), verts)
  vs = reinterpret(reshape, Float64, vs)
  modwrap = (x,a) -> mod(x-1, a) + 1
  a,c,d=acd
  nsides=length(verts)
  linear_polynomials=[[0,0],[1,0],[0,1]]
  H[1,1] = a
  for j=2:npolys
    for k=2:npolys
      for v=1:nsides
        vert = vs[:,v]
        next = vs[:,modwrap(v-1,nsides)]
        H[j,k] = H[j,k] + _integrate_scaled_monomials_on_triangle([c, vert, next], linear_polynomials[j], linear_polynomials[k], d)
      end
    end
  end
  H
end


function _generate_stabilising_term(Π∇, acd)
  B, D = Π∇
  projector = (B*D)\B
  nsides = size(B,2)
  (I(nsides) - D*projector)'*(I(nsides) - D*projector)
end

"""
For H matrix: Split polygons and integrate
"""
function _integrate_scaled_monomials_on_triangle(C, α, β, diameter)
  P1=C[1]
  P2=C[2]
  P3=C[3]

  centroid = P1
  M = [(1., Tuple(P1)...), (1., Tuple(P2)...), (1., Tuple(P3)...)]
  M = reinterpret(reshape, Float64, M)'
  area = 0.5*abs(det(M))

  Q = quadrature_rule(2)
  qw = Q[1]
  qx = Q[2]
  qy = Q[3]

  V = 0
  for q=1:length(qx)
    xhat = (P2[1]-P1[1])*qx[q] + (P3[1]-P1[1])*qy[q] + P1[1];
    yhat = (P2[2]-P1[2])*qx[q] + (P3[2]-P1[2])*qy[q] + P1[2];
    ma = ((xhat-centroid[1])/diameter)^(α[1])*((yhat-centroid[2])/diameter)^(α[2]);
    mb = ((xhat-centroid[1])/diameter)^(β[1])*((yhat-centroid[2])/diameter)^(β[2]);
    V = V + qw[q]*2*area*ma*mb
  end
  V
end
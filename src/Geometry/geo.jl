"""
Function to get the area, centroid, diameter of an element
"""

function geo(mesh::Triangulation, id::Int64)
  cell_verts = get_cell_coordinates(mesh)
  verts = cell_verts[id]
  if(length(verts) == 4)
    verts = view(verts,[1,2,4,3]) # To make it a closed quad
  end
  nsides = length(verts)
  vs = map(nc -> SVector(Tuple(nc)), verts)
  vs = reinterpret(reshape, Float64, vs)

  # Area
  area_comp = vs[1,:].*vs[2,vcat(2:end,1)] - vs[1, vcat(2:end,1)].*vs[2,:]
  area = 0.5*abs(sum(area_comp))

  # Centroid
  T₁ = transpose(repeat(area_comp, 1, 2))
  T₂ = (vs + vs[:, vcat(2:end,1)])
  T₁T₂ = (T₂.*T₁)/(6*area)
  centroid = Point(sum(T₁T₂, dims=2))

  # Diameter
  diameter=0
  for i=1:nsides-1
    for j=(i+1):nsides
      diameter = max(diameter, norm(vs[:,i] - vs[:,j], 2))
    end
  end

  area, centroid, diameter
end

function geo(mesh::Triangulation)
  lazy_map(x-> geo(mesh, x), 1:num_cells(mesh))
end

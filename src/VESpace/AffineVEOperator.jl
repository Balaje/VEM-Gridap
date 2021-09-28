import Gridap.FESpaces: SparseMatrixAssembler

"""
AffineVEOperator
"""

struct AffineVEOperator <: FEOperator
  trial::VESpace
  test::VESpace
  op::AffineOperator
end

function AffineVEOperator(trial::VESpace, test::VESpace, matrix::AbstractMatrix, vector::AbstractVector)
  @assert num_free_dofs(trial) == size(matrix,2) "Incompatible trial space and matrix"
  @assert num_free_dofs(test) == size(matrix,1) "Incompatible test space and matrix"
  op = AffineOperator(matrix,vector)
  AffineVEOperator(trial,test,op)
end

function AffineVEOperator(a::Function, l::Function, trial::VESpace, test::VESpace)
  mesh = Triangulation(trial.space.model)
  acd = geo(mesh)

  matcontribs = lazy_map(i -> _generate_mat_contribs(a, trial, acd, i), 1:num_cells(mesh))
  loadcontribs = lazy_map(i -> _generate_vec_contribs(l, trial, acd, i), 1:num_cells(mesh))

  ## Assemble the system
  σₖ = get_cell_dof_ids(trial.linear_fespace)
  matdata = ([matcontribs],[σₖ],[σₖ])
  vecdata = ([loadcontribs], [σₖ])
  assem = SparseMatrixAssembler(trial, test)
  matrix = assemble_matrix(assem, matdata)
  vector = assemble_vector(assem, vecdata)

  # Construct the Affine Operator
  op = AffineOperator(matrix, vector)
  AffineVEOperator(trial, test, op)
end

# Let us just return FEFunction
function Algebra.solve(op::AffineVEOperator)
  K = op.op.matrix
  f = op.op.vector
  u = K\f
  FEFunction(op.trial.linear_fespace, u)
end

function _generate_mat_contribs(a::Function, V::VESpace, acd, ind)
  ar,c,d = acd[ind]
  npolys = 3
  mesh = get_triangulation(V)
  stab_coeff = get_stab_coeff(V)

  H = zeros(npolys, npolys)
  linear_polynomials=[[0,0],[1,0],[0,1]]
  for i = 1:npolys
    polyᵢ = linear_polynomials[i]
    mᵢ(x) = ((x[1]- c[1])/d)^polyᵢ[1]*((x[2]- c[2])/d)^polyᵢ[2]
    cfᵢ = CellField(mᵢ, mesh)
    for j = 1:npolys
      polyⱼ = linear_polynomials[j]
      mⱼ(x) = ((x[1]- c[1])/d)^polyⱼ[1]*((x[2]- c[2])/d)^polyⱼ[2]
      cfⱼ = CellField(mⱼ, mesh)

      H[i,j] = a(cfᵢ, cfⱼ)[ind]
    end
  end

  B, D = get_Π∇(V)[ind]
  stability_term = get_stability(V)[ind]
  projector = (B*D)\B

  # Local matrix
  projector'*H*projector + sum(ar.^stab_coeff)*stability_term
end

function _generate_vec_contribs(l::Function, V::VESpace, acd, ind)
  ar, c, d = acd[ind]
  npolys = 3
  mesh = get_triangulation(V)

  Fⱼ = zeros(npolys, 1)
  linear_polynomials=[[0,0],[1,0],[0,1]]
  for i = 1:npolys
    polyᵢ = linear_polynomials[i]
    mᵢ(x) = ((x[1]- c[1])/d)^polyᵢ[1]*((x[2]- c[2])/d)^polyᵢ[2]
    cfᵢ = CellField(mᵢ, mesh)

    Fⱼ[i] = l(cfᵢ)[ind]
  end

  B, D = get_Π∇(V)[ind]
  projector = (B*D)\B

  # Local vector
  projector'*Fⱼ
end

function SparseMatrixAssembler(trial::VESpace, test::VESpace)
  SparseMatrixAssembler(trial.linear_fespace, test.linear_fespace)
end

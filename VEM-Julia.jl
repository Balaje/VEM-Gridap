using Gridap
using StaticArrays
using LinearAlgebra
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Algebra

include("src/Geometry/geo.jl")

include("src/VESpace/ConformingVESpace.jl")
include("src/VESpace/AffineVEOperator.jl")

include("src/Projectors/gradprojectors.jl")

include("src/Quadrature/quadrature_rule.jl")

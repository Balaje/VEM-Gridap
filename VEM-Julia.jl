using Gridap
using StaticArrays
using LinearAlgebra
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Algebra

include("src/Geometry/geo.jl")

include("src/VESpace/ConformingVESpace.jl")
include("src/VESpace/AffineVEOperator.jl")

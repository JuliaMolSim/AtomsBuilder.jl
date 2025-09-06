module AtomsBuilder

using AtomsBase 
using AtomsSystems
using LinearAlgebra: norm, I 
using Reexport
using StaticArrays: SMatrix, SVector
using StatsBase: sample, AbstractWeights
using Unitful

@reexport using StatsBase: AnalyticalWeights, FrequencyWeights, ProbabilityWeights, Weights

export bulk
export rattle! 
export randz!
export random_species!

const Mat3{T} = SMatrix{3, 3, T}
const Vec3{T} = SVector{3, T}

include("chemistry.jl")
include("utils.jl")
include("bulk.jl")
include("pubchem.jl")

end

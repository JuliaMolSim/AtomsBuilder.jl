module AtomsBuilder

import AtomsBase
using StaticArrays: SMatrix, SVector

const Mat3{T} = SMatrix{3, 3, T} 
const Vec3{T} = SVector{3, T}

include("chemistry.jl")

include("utils.jl")

include("bulk.jl")

end

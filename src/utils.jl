
using AtomsBase 
using AtomsBase: Atom, FlexibleSystem, Periodic
using Unitful: unit, ustrip 
using LinearAlgebra: norm 

export rattle!

"""
Helper function to convert construct a FlexibleSystem from 
   a list of positions, elements, cell matrix and pbc tuple 
"""
function _flexible_system(positions, elements, cell, pbc)
   Nat = length(positions)
   syms = Chemistry.chemical_symbol.(elements)
   atoms = [ Atom(; position = positions[i],
                    atomic_symbol = syms[i]) for i in 1:Nat ]
   bc =  [ (pbc[i] ? Periodic() : nothing) for i = 1:3 ]
   bb = [cell[i, :] for i = 1:3]
   return FlexibleSystem(atoms; 
                         bounding_box = bb, 
                         boundary_conditions = bc)
end

_set_position(x::Atom, 𝐫) = Atom(; position = 𝐫, 
                                   velocity = x.velocity,
                                   atomic_mass = x.atomic_mass,
                                   atomic_number = x.atomic_number, 
                                   atomic_symbol = x.atomic_symbol)

function _get_positions(at::FlexibleSystem)
   [ x.position for x in at.particles ]
end

function _set_positions(at::FlexibleSystem, 
                        X::AbstractVector{SVector{3, T}}) where {T}
   particles = [ _set_position(at.particles[i], X[i]) 
                 for i in 1:length(at) ] 
   return FlexibleSystem(particles, at.bounding_box, at.boundary_conditions)
end

function _get_atomic_numbers(at::FlexibleSystem)
   [ x.atomic_number for x in at.particles ]
end

function _set_atomic_numbers(at::FlexibleSystem, Z::AbstractVector{<: Integer})
   particles = [ Atom(Z[i], x.position, x.velocity)
                 for (i, x) in enumerate(at.particles) ]
   return FlexibleSystem(particles, at.bounding_box, at.boundary_conditions)
end


"""
```
repeat(at, n::NTuple{3})
repeat(at, n::Integer)
```

Takes a structure and repeats it n_j times
into the j-th cell-vector direction. For example,
```
at = repeat(bulk(:C), (3,2,4))
```
creates 3 x 2 x 4 unit cells of carbon.

The same can be achieved by `*`:
```
at = bulk(:) * (3, 2, 4)
```
"""
function Base.repeat(at::FlexibleSystem, n::NTuple{3})
   c1, c2, c3 = at.bounding_box

   particles = eltype(at.particles)[] 
   for a in CartesianIndices( (1:n[1], 1:n[2], 1:n[3]) )
      b = c1 * (a[1]-1) + c2 * (a[2]-1) + c3 * (a[3]-1)
      for i in 1:length(at)
         p_i = at.particles[i]
         p_new = _set_position(p_i, b + p_i.position)
         push!(particles, p_new)
      end
   end

   bb = [c1 * n[1], c2 * n[2], c3 * n[3]]
   return FlexibleSystem(particles, bb, at.boundary_conditions)
end

Base.repeat(at::FlexibleSystem, n::Integer) = repeat(at, (n,n,n))

import Base.*
*(at::FlexibleSystem, n) = repeat(at, n)
*(n, at::FlexibleSystem) = repeat(at, n)




"""
`rattle!(at, r::Float64) -> at`

Randomly perturbs the atom positions within a ball of radius `r`. The perturbation 
is uniform in angular component, and uniform in radial component. (Note this is 
not the same as choosing them uniform in cartesian coordinates!). 
"""
function rattle!(at::FlexibleSystem, r::AbstractFloat)
   for i = 1:length(at.particles)
      p = at.particles[i]
      𝐫i = p.position 
      T = typeof(ustrip(𝐫i[1]))
      ui = randn(Vec3{T})
      p_new = _set_position(p, 𝐫i + r * ui / norm(ui) * unit(𝐫i[1]))
      at.particles[i] = p_new
   end
   return at
end



# union(at1::Atoms, at2::Atoms) =
#    Atoms( X = union(at1.X, at2.X),
#           P = union(at1.P, at2.P),
#           M = union(at1.M, at2.M),
#           Z = union(at1.Z, at2.Z),
#           cell = cell(at1),
#           pbc = pbc(at1) )

# append(at::Atoms, X::AbstractVector{<:JVec}) =
#    Atoms( X = union(at.X, X),
#           P = union(at.P, zeros(JVecF, length(X))),
#           M = union(at.M, zeros(length(X))),
#           Z = union(at.Z, zeros(Int, length(X))),
#           cell = cell(at),
#           pbc = pbc(at) )


# """
# `deleteat!(at::Atoms, n) -> at`:

# returns the same atoms object `at`, but with the atom(s) specified by `n`
# removed.
# """
# function Base.deleteat!(at::Atoms, n)
#    deleteat!(at.X, n)
#    deleteat!(at.P, n)
#    deleteat!(at.M, n)
#    deleteat!(at.Z, n)
#    update_data!(at, Inf)
#    JuLIP.reset_clamp!(at)
#    return at
# end



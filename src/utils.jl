


"""
Helper function to convert construct a FlexibleSystem from 
   a list of positions, elements, cell matrix and pbc tuple 
"""
function _flexible_system(positions, elements, cell, pbc)
   Nat = length(positions)
   syms = Chemistry.chemical_symbol.(elements)
   atoms = [ Atom(syms[i], positions[i]) for i in 1:Nat ]
   return FlexibleSystem(atoms;
                         cell_vectors=tuple([cell[i, :] for i = 1:3]...),
                         periodicity=pbc)
end

_set_position(x::Atom, 𝐫) = Atom(atomic_number(x), 𝐫; 
                                   velocity = velocity(x),
                                   mass = mass(x))

_set_element(x::Atom, Z) = Atom(Z, position(x); 
                                velocity = velocity(x),
                                mass = mass(x))
 
function set_positions(at::FlexibleSystem, 
                       X::AbstractVector{SVector{3, T}}) where {T}
   @assert length(X) == length(at)                       
   particles = [ _set_position(at.particles[i], X[i])
                 for i in 1:length(at) ] 
   return FlexibleSystem(particles, 
                         cell_vectors = cell_vectors(at), 
                         periodicity = periodicity(at))
end


function set_elements(at::FlexibleSystem, Z::AbstractVector)
   @assert length(Z) == length(at)                       
   particles = [ Atom(Z[i], position(x), velocity(x))
                 for (i, x) in enumerate(at.particles) ]
   return FlexibleSystem(particles, 
                         cell_vectors = cell_vectors(at), 
                         periodicity = periodicity(at))
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
   c1, c2, c3 = cell_vectors(at)

   particles = eltype(at.particles)[] 
   for a in CartesianIndices( (1:n[1], 1:n[2], 1:n[3]) )
      b = c1 * (a[1]-1) + c2 * (a[2]-1) + c3 * (a[3]-1)
      for i in 1:length(at)
         p_i = at.particles[i]
         p_new = _set_position(p_i, b + p_i.position)
         push!(particles, p_new)
      end
   end

   bb = (c1 * n[1], c2 * n[2], c3 * n[3])
   return FlexibleSystem(particles; cell_vectors = bb, 
                           periodicity = periodicity(at))
end

Base.repeat(at::FlexibleSystem, n::Integer) = repeat(at, (n,n,n))

import Base.*
*(at::FlexibleSystem, n) = repeat(at, n)
*(n, at::FlexibleSystem) = repeat(at, n)


@deprecate *(sys::CellSystem, n) repeat(sys, n)
@deprecate *(n, sys::CellSystem) repeat(sys, n)

@deprecate rattle!(sys::CellSystem, r::Union{AbstractFloat, Quantity}) rattle_positions!(sys,r)

"""
```
rattle!(sys, r::Union{AbstractFloat, Quantity}) -> at
```

Randomly perturbs the atom positions within a ball of radius `r`. The perturbation 
is uniform in angular component, and uniform in radial component. (Note this is 
not the same as choosing them uniform in cartesian coordinates!). 

If `r` is unitless, then the unit of the system is applied. 
"""
function rattle!(at::FlexibleSystem, r::Quantity)
   (length(at.particles) > 0) || return at
   at_unit = unit(position(at, 1)[1])
   r = uconvert(at_unit, r)
   for i = 1:length(at.particles)
      p = at.particles[i]
      𝐫ᵢ = p.position 
      T = typeof(ustrip(𝐫ᵢ[1]))
      ui = randn(Vec3{T})
      p_new = _set_position(p, 𝐫ᵢ + rand(T) * r * ui / norm(ui))
      at.particles[i] = p_new
   end
   return at
end

rattle!(sys::FlexibleSystem, r::AbstractFloat) = 
      rattle!(sys, r * unit(position(sys, 1)[1]))


"""
   randz!(sys::FlexibleSystem, zlist) -> sys

Randomly assigns elements to the atoms in the system `sys` according to the
probabilities given in `zlist`.
`zlist` is an iterable over pairs of the form `id => p` where `id`
is an atom id (e.g. atomic number or chemical symbol) and `p` 
a probability. E.g., 

```julia
sys = bulk(:Ti, cubic=true) * 3
sys = randz!(sys, [ :Ti => 0.2, :O => 0.8 ])
```

This function was developed mostly for generating testing 
systems. It may not be suitable for generating random alloys. 
PRs to improve it are welcome. 
"""
function randz!(sys::AbstractSystem, zlist)
   species = [ x[1] for x in zlist ]
   weights = [ x[2] for x in zlist ]
   return random_species!(sys, species, Weights(weights))
end

"""
   `random_species!(sys::AbstractSystem, species::AbstractVector, weights)`

   `random_species!(sys::AbstractSystem, species::AbstractVector{ChemicalSpecies}, weights)`

Randomly assigns elements to the atoms in the system `sys` according to the
probabilities given in `weights`.

# Arguments
- `sys::AbstractSystem`: the system to modify
- `species::AbstractVector`: a vector of species (either atomic numbers or chemical symbols)
- `weights`: a `AbstractWeights` object from `StatsBase` giving the probabilities for each species

# Example
```julia
using AtomsBuilder

sys = repeat(bulk(:Ti), 3)

spc = [:Ti, :Al]  # species to sample from

# Create weights for each species
w = AnalyticalWeights([0.5, 0.5])
w = FrequencyWeights([1, 3])
w = ProbabilityWeights([0.2, 0.8])

random_species!( sys, spc, w )
```
"""
function random_species!(sys::AbstractSystem, species::AbstractVector, weights)
   spc = ChemicalSpecies.(species)
   return random_species!(sys, spc, weights)
end

function random_species!(sys::AbstractSystem, species::AbstractVector{ChemicalSpecies}, weights::AbstractWeights)
   for i in 1:length(sys)
      z = sample(species, weights)
      AtomsBase.set_species!(sys, i, z)
   end
   return sys
end


"""
```
union(sys1::FlexibleSystem, sys2::FlexibleSystem)
```
takes the union of two particle systems provided their cells are identical. 
"""
function union(sys1::FlexibleSystem, sys2::FlexibleSystem) 
   @assert periodicity(sys1) == periodicity(sys2)
   @assert cell_vectors(sys1) == cell_vectors(sys2)
   return FlexibleSystem(union(sys1.particles, sys2.particles),  
                        cell_vectors = cell_vectors(at), 
                        periodicit = periodicity(at) )
end

@deprecate union(sys1::CellSystem, sys2::CellSystem) union(sys1, sys2)

"""
`deleteat!(sys::FlexibleSystem, n) -> sys`:

returns the same `FlexibleSystem`` object `sys`, but with the atom(s) specified by `n`
removed.
"""
function Base.deleteat!(sys::FlexibleSystem, n)
   deleteat!(sys.particles, n)
   return sys
end 

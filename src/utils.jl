
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
w = FrequencyWeights([1, 3])
w = ProbabilityWeights([0.2, 0.8])
w = Weights([0.2, 0.8]) # generic weights

random_species!( sys, spc, w )
```
"""
function random_species!(sys::AbstractSystem, species::AbstractVector, weights)
   spc = ChemicalSpecies.(species)
   return random_species!(sys, spc, weights)
end

function random_species!(
      sys::AbstractSystem,
      species::AbstractVector{ChemicalSpecies},
      weights::StatsBase.AbstractWeights
   )
   for i in 1:length(sys)
      z = sample(species, weights)
      AtomsBase.set_species!(sys, i, z)
   end
   return sys
end

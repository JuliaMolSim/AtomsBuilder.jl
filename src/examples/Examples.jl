"""
    AtomsBuilder.Examples

Staging area for builders of example / fixture systems that are
useful to downstream packages for testing, benchmarking, and tuning.

Currently exposes:

- [`rocksalt`](@ref) — 1:1 rock-salt supercell of two species (NaCl
  structure), with charges as a per-atom property.
- [`nacl`](@ref)     — convenience alias for `rocksalt(:Na, :Cl, n)`.
- [`tip3p_water`](@ref) — cubic box of randomly-placed, randomly-oriented
  TIP3P water molecules.

These builders are kept in a submodule (rather than top-level) because
the API may evolve as we learn what consumers want; pin a version if
you depend on them.

Usage:

```julia
using AtomsBuilder.Examples
sys = nacl(2)                     # 64-ion NaCl supercell
sys = tip3p_water(12.0u"Å")       # ~48-molecule water box
```
"""
module Examples

using AtomsBase: Atom, FlexibleSystem
using Random: AbstractRNG, default_rng, randn, rand
using Unitful
using UnitfulAtomic               # provides u"e_au" (elementary charge)
using LinearAlgebra: norm
using StaticArrays: SVector, SMatrix

using ..AtomsBuilder: Vec3, Mat3, _convert_pbc

export rocksalt, nacl, tip3p_water

include("rocksalt.jl")
include("poisson.jl")
include("water.jl")

end  # module Examples

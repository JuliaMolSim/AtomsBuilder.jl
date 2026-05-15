# AtomsBuilder

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/AtomsBuilder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/AtomsBuilder.jl/dev/)-->
[![Build Status](https://github.com/JuliaMolSim/AtomsBuilder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaMolSim/AtomsBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain) 

This package provides utilities to build atomic structures. At the moment the functionality is limited - see examples below. The intention is that over time this package becomes the de facto standard for generating structures in the JuliaMolSim ecosystem. Contributions are very welcome.


## Preliminary Documentation 

Currently there are just two exported functions to build materials: 
* `bulk`
* `rattle!`
In addition we overload 
* `repeat` (with alias `*`)

```julia
using AtomsBuilder 

# generate a diamond cubic bulk Si unit cell 
at1 = bulk(:Si)
@show length(at1)

# generate a minimal cubic Si cell (diamond cubic)
at2 = bulk(:Si, cubic=true)
@show length(at2)

# repeat the cell 3 times in each coordinate direction
at3 = at2 * 3
@show length(at3)

# repeat the unit cell in only one direction
at4 = at2 * (3, 1, 1)
@show length(at3)

# create a bulk supercell and then rattle the atoms 
at5 = rattle!( bulk(:Si, cubic=true) * 3 )
```

See `?bulk` and `?rattle!` for more information. 

## `AtomsBuilder.Examples`

The `Examples` submodule ships fixture systems that are useful for
testing, benchmarking, and tuning downstream packages. The API may
evolve — pin a version if you depend on it.

* `rocksalt(species_a, species_b, n)` — 1:1 rock-salt supercell (NaCl
  structure) with charges `±1 e_au` attached as a per-atom property.
* `nacl(n)` — convenience alias for `rocksalt(:Na, :Cl, n)`.
* `tip3p_water(box)` — cubic box of randomly-placed, randomly-oriented
  TIP3P water molecules (Bridson Poisson-disk for O placement; charges
  `q_O = -0.834 e_au`, `q_H = +0.417 e_au`).

```julia
using AtomsBuilder.Examples, Unitful, UnitfulAtomic

# 64-ion NaCl supercell, perfect lattice
sys = nacl(2)

# 216-ion NaCl supercell, σ = 0.1 Å Gaussian displacements (~300 K thermal)
sys = nacl(3; σ = 0.1u"Å", rng = Random.MersenneTwister(0))

# ~58-molecule water box, default density (1 g/cm³), default d_min = 2.7 Å
sys = tip3p_water(12.0u"Å")
```

Charges are stored on each `Atom` as a per-atom property under the
keyword `:charge` by default (extracted via `atom.data[:charge]`).
Pass `charge_label = :q` (or similar) to use a different key.

See `?rocksalt`, `?nacl`, `?tip3p_water` for full options.

## PubChem Interface

PubChem interface allows you to download structures from [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
This is done with `load_from_pubchem` function. You can only load isolated molecules currently.

```julia
using AtomsBuilder

# using trivial name
load_from_pubchem( "water" )

# using CID
load_from_pubchem( 887 )

# using SMILES
load_from_pubchem( smiles="CC(=O)C" )

# using CAS number
load_from_pubchem( "64-17-5" )
```

## Contributions 

The the package started as a copy-paste of a subset of functionality from an older package that is no longer developed. Contributions to expand the capabilities, improve the implementation, or entirely replace it are very welcome. There are almost certainly more general and more elegant implementations of structure building available than what we currently. Some packages that contain overlapping functionalities that could replace or add to `AtomsBuilder.jl` include
* [`Electrum.jl`](https://github.com/brainandforce/Electrum.jl)
* [`AtomsToolbox.jl`](https://github.com/rashidrafeek/AtomsToolbox.jl)
* [`SimpleCrystals.jl`](https://github.com/ejmeitz/SimpleCrystals.jl)
* [`Packmol.jl`](https://github.com/m3g/Packmol.jl)
* [`ASEconvert.jl`](https://github.com/mfherbst/ASEconvert.jl)


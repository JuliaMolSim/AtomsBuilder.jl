# AtomsBuilder

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/AtomsBuilder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/AtomsBuilder.jl/dev/)-->
[![Build Status](https://github.com/JuliaMolSim/AtomsBuilder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaMolSim/AtomsBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain) 

This package provides utilities to build atomic structures. At the moment the functionality is very limited - see examples below. 


## Preliminary Documentation 

Currently there are just two exported functions: 
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

## Contributions 

The current version of the package is essentially a copy-paste of a subset of functionality from an older package that is no longer developed. Contributions to expand the capabilities, improve the implementation, or entirely replace it are very welcome. Some packages that contain overlapping functionalities that could replace or add to `AtomsBuilder.jl` include
* [`Electrum.jl`](https://github.com/brainandforce/Electrum.jl)
* [`AtomsToolbox.jl`](https://github.com/rashidrafeek/AtomsToolbox.jl)

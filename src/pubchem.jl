using AtomsBase: FastSystem
using AtomsIO: load_system
using PubChemCrawler: get_cid, get_for_cids 

export load_from_pubchem

"""
    load_from_pubchem(cid::Int)
    load_from_pubchem(name::AbstractString)
    load_from_pubchem(;smiles::AbstractString)

Load a molecule from [PubChem](https://pubchem.ncbi.nlm.nih.gov/).

You can use [PubChemCrawler](https://github.com/JuliaHealth/PubChemCrawler.jl)
to search for `cid`.

# Example
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
"""
function load_from_pubchem(cid::Int)
    tmp_name = tempname() * ".sdf"
    open(tmp_name, "w") do io
        write(io, get_for_cids(cid, output="SDF", record_type="3d"))
    end
    # parsing system does not work properly so we drop global features
    sys = FastSystem( load_system(tmp_name) )
    rm(tmp_name)
    return sys
end

function load_from_pubchem(name::AbstractString)
    cid = get_cid(name=name)
    return load_from_pubchem(cid)
end

function load_from_pubchem(;smiles::AbstractString)
    cid = get_cid(smiles=smiles)
    return load_from_pubchem(cid)
end
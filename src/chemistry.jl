module Chemistry

using JSON 

# --------------------------------------------------------------------
# Load data and prepare it a bit ... 

const ase_data_path = joinpath(@__DIR__(), "..", "data", "asedata.json")

const ase_data = JSON.parsefile(ase_data_path)

const _rnn = [Float64(d) for d in ase_data["rnn"]]
const _masses = [Float64(m) for m in ase_data["masses"]]
const _refstates = Dict{String, Any}[ d == nothing ? Dict{String, Any}() : d
                                      for d in ase_data["refstates"]]

const _symbols = Dict{Int, Symbol}()
const _numbers = Dict{Symbol, Int}()
for (n, sym) in enumerate(Symbol.(ase_data["symbols"]))
   _numbers[sym] = Int(n - 1)
   _symbols[n-1] = sym
end

# --------------------------------------------------------------------
#  some convenient accessor and conversion functions 


atomic_number(sym::Symbol) = _numbers[sym]

chemical_symbol(z::Integer) = _symbols[z]

atomic_mass(z::Integer) = _masses[z+1]
atomic_mass(sym::Symbol) = atomic_mass(atomic_number(sym))

rnn(z::Integer) = _rnn[z+1]
rnn(sym::Symbol) = rnn(atomic_number(sym))

symmetry(z::Integer) = Symbol(_refstates[z+1]["symmetry"])
symmetry(sym::Symbol) = symmetry(atomic_number(sym))

lattice_constant(z::Integer) = Float64( _refstates[z+1]["a"] )
lattice_constant(sym::Symbol) = lattice_constant(atomic_number(sym))

refstate(z::Integer) = _refstates[z+1]
refstate(sym::Symbol) = refstate(atomic_number(sym))

end
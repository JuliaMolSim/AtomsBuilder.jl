export bulk

using Unitful: unit
using LinearAlgebra: I

const _unit_cells = Dict(     # (positions, cell matrix, factor of a)
   :fcc => ( [ [0 0 0], ], [0 1 1; 1 0 1; 1 1 0],  0.5),
   :bcc => ( [ [0.0,0.0,0.0], ],
             [-1 1 1; 1 -1 1; 1 1 -1], 0.5),
   :diamond => ( [ [0.0, 0.0, 0.0], [0.5, 0.5, 0.5] ],
                 [0 1 1; 1 0 1; 1 1 0], 0.5)
)

const _cubic_cells = Dict(   # (positions, factor of a)
   :fcc     => ( [ [0 0 0], [0 1 1], [1 0 1], [1 1 0] ], 0.5 ),
   :bcc     => ( [ [0 0 0], [1 1 1] ], 0.5 ),
   :diamond => ( [ [0 0 0], [1 1 1], [0 2 2], [1 3 3], [2 0 2],
                   [3 1 3], [2 2 0], [3 3 1] ], 0.25 )
)

function _simple_bulk(sym::Symbol, cubic::Bool; a=nothing, T=Float64)
   if cubic
      X, scal = _cubic_cells[Chemistry.symmetry(sym)]
      C = Matrix(1.0I, 3, 3) / scal
   else
      X, C, scal = _unit_cells[Chemistry.symmetry(sym)]
   end
   a = @something a Chemistry.lattice_constant(sym)
   TU = typeof( one(T) * unit(a) )
   [ Vec3{TU}(x * a * scal)  for x in X ], Mat3{TU}(C * a * scal)
end


function _bulk_hcp(sym::Symbol; a=nothing, c=nothing, T=Float64)
   D = Chemistry.refstate(sym)
   a = @something a   D["a"]
   c = @something c a*D["c/a"]
   [a * Vec3{T}(0.0,         0.0,       0.0),
    a * Vec3{T}(0.0, 1 / sqrt(3), c / a / 2) ],
    a * Mat3{T}([1,      -1/2,   0,
                 0, sqrt(3)/2,   0,
                 0,         0, c/a])
end


_convert_pbc(pbc::NTuple{3, Bool}) = pbc
_convert_pbc(pbc::Bool) = (pbc, pbc, pbc)
_convert_pbc(pbc::AbstractVector) = tuple(pbc...)

"""
`bulk(sym)` : generates a `FlexibleSystem` unit cell for a bulk
crystal structure. If `sym` is a chemical symbol then the phase and
lattice constant are taken from a database that is consistent with ASE.
Optional alternative values can be chosen via the kwargs
`a`, `b` or `c` to specify the lattice constants.
"""
function bulk(sym::Symbol; cubic=false, pbc=(true, true, true),
                           a=nothing, c=nothing, b=nothing, T=Float64)
                           # , x=nothing, y=nothing, z=nothing)
    symm = try
        Chemistry.symmetry(sym)
    catch e
        if e isa KeyError
            throw(ArgumentError("No symmetry information known for element $sym"))
        else
            rethrow()
        end
    end
    if symm in (:fcc, :bcc, :diamond)
        X, C = _simple_bulk(sym, cubic; a, T)
    elseif symm == :hcp
        cubic && throw(ArgumentError("cubic=true cannot be selected for $symm lattices"))
        X, C = _bulk_hcp(sym; a, c, T)
    else
        throw(ArgumentError("Currently bulk not implemented for symmetry $symm"))
    end
    #Z = Chemistry.atomic_number(sym)
    nat = length(X)
    # at = _flexible_system( X, fill(Z, nat), C, _convert_pbc(pbc))
    spc = AtomsBase.ChemicalSpecies(sym)
    at = generic_system( fill(spc, nat), X; cell_vectors=collect(eachrow(C)), periodicity=pbc)
    # I don't remember what this does so I'm commenting it out until I understand 
    # it again ...                         
    # (x !== nothing || y !== nothing || z !== nothing) && rotate!(at, x=x, y=y, z=z)
    #return at
    return at
end

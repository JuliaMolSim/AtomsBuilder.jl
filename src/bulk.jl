export bulk

using Unitful: unit 
using LinearAlgebra: I

const _unit_cells = Dict(     # (positions, cell matrix, factor of a)
   :fcc => ( [ [0 0 0], ],
             [0 1 1; 1 0 1; 1 1 0],  0.5),
   :bcc => ( [ [0.0,0.0,0.0], ],
             [-1 1 1; 1 -1 1; 1 1 -1], 0.5),
   :diamond => ( [ [0.0, 0.0, 0.0], [0.5, 0.5, 0.5] ],
                 [0 1 1; 1 0 1; 1 1 0], 0.5)
)

const _cubic_cells = Dict(   # (positions, factor of a)
   :fcc => ( [ [0 0 0], [0 1 1], [1 0 1], [1 1 0] ], 0.5 ),
   :bcc => ( [ [0 0 0], [1 1 1] ], 0.5 ),
   :diamond => ( [ [0 0 0], [1 1 1], [0 2 2], [1 3 3], [2 0 2],
                   [3 1 3], [2 2 0], [3 3 1] ], 0.25 )
)

const _simple_structures = [:fcc, :bcc, :diamond]

function _simple_bulk(sym::Symbol, cubic::Bool; a=nothing, T = Float64)
   if cubic
      X, scal = _cubic_cells[Chemistry.symmetry(sym)]
      C = Matrix(1.0I, 3, 3) / scal
   else
      X, C, scal = _unit_cells[Chemistry.symmetry(sym)]
   end
   a === nothing && (a = Chemistry.lattice_constant(sym))
   TU = typeof( one(T) * unit(a) )
   return [ Vec3{TU}(x * a * scal)  for x in X ], Mat3{TU}(C * a * scal)
end


function _bulk_hcp(sym::Symbol; a=nothing, c=nothing, T=Float64)
   D = Chemistry.refstate(sym)
   a === nothing && (a = D["a"])
   c === nothing && (c = a * D["c/a"])
   return [ a * Vec3{T}(0.0, 0.0, 0.0), 
            a * Vec3{T}(0.0, 1 / sqrt(3), c / a / 2) ],
         a * Mat3{T}( [1.0, -1/2,  0.0, 0.0, sqrt(3)/2, 0.0, 0.0, 0.0, c/a] )
end


_convert_pbc(pbc::NTuple{3, Bool}) = pbc
_convert_pbc(pbc::Bool) = (pbc, pbc, pbc)

"""
`bulk(sym)` : generates a `FlexibleSystem` unit cell for a bulk 
crystal structure. If `sym` is a chemical symbol then the phase and  
lattice constant are taken from a database that is consistent with ASE. 
If `sym` is one of `[:fcc, :bcc, :diamond, :hcp]` then one needs to 
specify the kwargs `a` or `c` to determine the lattice constants. 
"""
function bulk(sym::Symbol; T=Float64, cubic = false, pbc = (true,true,true), 
                           a=nothing, c=nothing) # , x=nothing, y=nothing, z=nothing)
   symm = Chemistry.symmetry(sym)
   if symm in _simple_structures
      X, C = _simple_bulk(sym, cubic; a=a)
   elseif symm == :hcp
      X, C = _bulk_hcp(sym; a=a, c=c)  # cubic parameter is irrelevant for hcp
   end
   m = Chemistry.atomic_mass(sym)
   Z = Chemistry.atomic_number(sym)
   nat = length(X)
   at = _flexible_system( X, fill(Z, nat), C, _convert_pbc(pbc))
   # I don't remember what this does so I'm commenting it out until I understand 
   # it again ...                         
   # (x !== nothing || y !== nothing || z !== nothing) && rotate!(at, x=x, y=y, z=z)
   return at
end

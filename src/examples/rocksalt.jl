# --- rock-salt supercell -----------------------------------------------
# Reference structure for 1:1 ionic crystals (NaCl, LiF, KCl, MgO, ...).
# The conventional cell hosts 4 species-A + 4 species-B at the standard
# fcc / interpenetrating-fcc fractional sites.

const _rocksalt_a_frac = (SVector(0.0, 0.0, 0.0), SVector(0.5, 0.5, 0.0),
                          SVector(0.5, 0.0, 0.5), SVector(0.0, 0.5, 0.5))
const _rocksalt_b_frac = (SVector(0.5, 0.0, 0.0), SVector(0.0, 0.5, 0.0),
                          SVector(0.0, 0.0, 0.5), SVector(0.5, 0.5, 0.5))

# Build an Atom with `charge_label => q` attached as data.
# Splatting a single Pair into kwargs lets `charge_label` be a runtime Symbol.
@inline _charged_atom(sym, pos, q, charge_label) =
   Atom(sym, pos; (charge_label => q,)...)

"""
    rocksalt(species_a, species_b, n_super;
             a_lat            = 4.0u"Å",
             σ                = 0.0u"Å",
             pbc              = (true, true, true),
             charge_label::Symbol = :charge,
             rng::AbstractRNG = default_rng())

Rock-salt (NaCl-structure) supercell of `species_a` and `species_b` in
a 1:1 ratio. `n_super` is the number of conventional unit cells per
axis, giving `8 · n_super³` ions in total.

Per-atom charges (`+1 e_au` for `species_a`, `-1 e_au` for `species_b`)
are attached as a per-atom property with key `charge_label` (default
`:charge`). `a_lat` is the conventional-cell lattice constant. If
`σ > 0`, every ion gets an independent isotropic Gaussian Cartesian
displacement of standard deviation `σ` — `σ ≈ 0.1 Å` roughly mimics a
300 K thermal sample of NaCl.

Bare-`Real` values for `a_lat` and `σ` are taken to be in Å.

```julia
sys = rocksalt(:Na, :Cl, 3)                        # 216-ion NaCl
sys = rocksalt(:Li, :F,  2; a_lat = 4.02u"Å")      # 64-ion LiF
sys = rocksalt(:Na, :Cl, 2; σ = 0.1u"Å",
                            rng = MersenneTwister(0))   # thermal
```
"""
function rocksalt(species_a::Symbol, species_b::Symbol, n_super::Integer;
                  a_lat                = 4.0u"Å",
                  σ                    = 0.0u"Å",
                  pbc                  = (true, true, true),
                  charge_label::Symbol = :charge,
                  rng::AbstractRNG     = default_rng())
   a_lat = a_lat isa Real ? a_lat * u"Å" : a_lat
   σ     = σ     isa Real ? σ     * u"Å" : σ

   TU       = typeof(a_lat)
   q_a      = +1.0u"e_au"
   q_b      = -1.0u"e_au"
   n_atoms  = 8 * n_super^3
   atoms    = Vector{Atom}(); sizehint!(atoms, n_atoms)
   thermal  = !iszero(σ)

   @inbounds for i in 0:n_super-1, j in 0:n_super-1, k in 0:n_super-1
      offset = SVector(Float64(i), Float64(j), Float64(k)) * a_lat
      for p in _rocksalt_a_frac
         r0 = p * a_lat + offset
         δ  = thermal ? σ * SVector(randn(rng), randn(rng), randn(rng)) :
                         zero(Vec3{TU})
         push!(atoms, _charged_atom(species_a, r0 + δ, q_a, charge_label))
      end
      for p in _rocksalt_b_frac
         r0 = p * a_lat + offset
         δ  = thermal ? σ * SVector(randn(rng), randn(rng), randn(rng)) :
                         zero(Vec3{TU})
         push!(atoms, _charged_atom(species_b, r0 + δ, q_b, charge_label))
      end
   end

   box  = n_super * a_lat
   cell = Mat3{TU}(box * one(SMatrix{3, 3, Float64}))
   return FlexibleSystem(atoms;
                         cell_vectors = tuple([cell[i, :] for i = 1:3]...),
                         periodicity  = _convert_pbc(pbc))
end

"""
    nacl(n_super; kwargs...)

Convenience alias for `rocksalt(:Na, :Cl, n_super; kwargs...)`. Defaults
match standard rock-salt NaCl (`a_lat = 4.0u"Å"` ≈ ASE reference).

```julia
sys = nacl(2)              # 64-ion NaCl supercell, perfect lattice
sys = nacl(3; σ = 0.1)     # 216-ion, 300 K thermal sample
```
"""
nacl(n_super::Integer; kwargs...) = rocksalt(:Na, :Cl, n_super; kwargs...)

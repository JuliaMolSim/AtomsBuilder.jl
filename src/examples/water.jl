# --- TIP3P-like water box ----------------------------------------------
# TIP3P (Jorgensen et al., J. Chem. Phys. 79, 926 (1983)): rigid 3-site
# water model. q_O = -0.834 e_au, q_H = +0.417 e_au, r_OH = 0.9572 Å,
# ∠HOH = 104.52°.

const _TIP3P_R_OH_VAL = 0.9572                         # Å (bare for math)
const _TIP3P_THETA    = 104.52 * π / 180
const _TIP3P_Q_O      = -0.834u"e_au"
const _TIP3P_Q_H      = +0.417u"e_au"

# H positions in the molecule's local frame (O at origin, C₂ axis ‖ ẑ).
# Stored as unitless Vec3{Float64} (Å implicit); the rotation matrix R
# is then unitless 3×3 and `R * H_local` is unitless 3-vector; units are
# applied when constructing the final Atom positions.
const _TIP3P_H1_LOCAL = SVector{3, Float64}(
   _TIP3P_R_OH_VAL * sin(_TIP3P_THETA / 2), 0.0,
   _TIP3P_R_OH_VAL * cos(_TIP3P_THETA / 2),
)
const _TIP3P_H2_LOCAL = SVector{3, Float64}(
  -_TIP3P_R_OH_VAL * sin(_TIP3P_THETA / 2), 0.0,
   _TIP3P_R_OH_VAL * cos(_TIP3P_THETA / 2),
)

"""
    tip3p_water(box;
                ρ                = 0.0334u"Å^-3",
                d_min            = 2.7u"Å",
                pbc              = (true, true, true),
                charge_label::Symbol = :charge,
                rng::AbstractRNG = default_rng())

Cubic-box TIP3P liquid-water configuration. Side length `box`, target
number density `ρ` ⇒ `n_mol = round(ρ · box³)` molecules. Oxygens are
placed by Bridson Poisson-disk sampling (min O–O distance `d_min`,
minimum-image PBC). Each molecule gets a uniformly-random SO(3)
orientation. Atom order per molecule is `[O, H, H]`.

Per-atom charges (`q_O = -0.834 e_au`, `q_H = +0.417 e_au`) are
attached as a per-atom property with key `charge_label` (default
`:charge`).

Bare-`Real` values for `box` and `d_min` are taken to be in Å; bare
`ρ` is in Å⁻³.

Standard liquid-water values at 300 K: `ρ = 0.0334 mol/Å³` ≈ 1 g/cm³,
`d_min ≈ 2.7 Å` (just below the first O–O peak in the RDF).

```julia
sys = tip3p_water(12.0u"Å")                           # ~58 molecules
sys = tip3p_water(15.0; rng = MersenneTwister(0))     # pinned
```
"""
function tip3p_water(box;
                     ρ                    = 0.0334u"Å^-3",
                     d_min                = 2.7u"Å",
                     pbc                  = (true, true, true),
                     charge_label::Symbol = :charge,
                     rng::AbstractRNG     = default_rng())
   box   = box   isa Real ? box   * u"Å"    : box
   ρ     = ρ     isa Real ? ρ     * u"Å^-3" : ρ
   d_min = d_min isa Real ? d_min * u"Å"    : d_min

   box_val   = ustrip(u"Å",    box)
   d_min_val = ustrip(u"Å",    d_min)
   ρ_val     = ustrip(u"Å^-3", ρ)

   n_mol     = round(Int, ρ_val * box_val^3)
   O_centres = _poisson_disk_oxygens(box_val, n_mol, d_min_val, rng)

   atoms = Vector{Atom}(); sizehint!(atoms, 3 * n_mol)
   for O in O_centres
      R  = _random_rotation_matrix(rng)
      H1 = O + R * _TIP3P_H1_LOCAL
      H2 = O + R * _TIP3P_H2_LOCAL
      push!(atoms, _charged_atom(:O, O  * u"Å", _TIP3P_Q_O, charge_label))
      push!(atoms, _charged_atom(:H, H1 * u"Å", _TIP3P_Q_H, charge_label))
      push!(atoms, _charged_atom(:H, H2 * u"Å", _TIP3P_Q_H, charge_label))
   end

   TU   = typeof(box)
   cell = Mat3{TU}(box * one(SMatrix{3, 3, Float64}))
   return FlexibleSystem(atoms;
                         cell_vectors = tuple([cell[i, :] for i = 1:3]...),
                         periodicity  = _convert_pbc(pbc))
end

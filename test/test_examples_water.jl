
using AtomsBuilder
using AtomsBuilder.Examples
using Test, AtomsBase, Unitful, UnitfulAtomic, Random
using LinearAlgebra: norm, ⋅

##

@testset "tip3p_water — size, density, structure" begin
   box  = 12.0u"Å"
   ρ    = 0.0334u"Å^-3"
   sys  = tip3p_water(box; ρ = ρ, rng = MersenneTwister(0))

   n_mol = round(Int, ustrip(u"Å^-3", ρ) * ustrip(u"Å", box)^3)
   @test length(sys) == 3 * n_mol
   @test length(sys) % 3 == 0
   @test periodicity(sys) == (true, true, true)

   # Cubic cell of side `box`.
   cv = cell_vectors(sys)
   for (k, v) in enumerate(cv)
      @test ustrip(u"Å", v[k]) ≈ ustrip(u"Å", box)
   end

   # Atom order is [O, H, H, O, H, H, ...].
   zs = atomic_number(sys, :)
   z_O = AtomsBuilder.Chemistry.atomic_number(:O)
   z_H = AtomsBuilder.Chemistry.atomic_number(:H)
   @test zs[1:3:end] == fill(z_O, n_mol)
   @test zs[2:3:end] == fill(z_H, n_mol)
   @test zs[3:3:end] == fill(z_H, n_mol)
end

@testset "tip3p_water — charges (TIP3P values, neutral molecules)" begin
   sys = tip3p_water(12.0u"Å"; rng = MersenneTwister(0))
   charges = [a.data[:charge] for a in sys]
   @test charges[1:3:end] == fill(-0.834u"e_au", length(sys) ÷ 3)
   @test charges[2:3:end] == fill(+0.417u"e_au", length(sys) ÷ 3)
   @test charges[3:3:end] == fill(+0.417u"e_au", length(sys) ÷ 3)
   # Each molecule is neutral to high precision; total exactly zero by symmetry.
   total = sum(charges)
   @test abs(ustrip(u"e_au", total)) < 1e-12
end

@testset "tip3p_water — O–H bond length and H–O–H angle" begin
   sys = tip3p_water(12.0u"Å"; rng = MersenneTwister(0))
   X   = position(sys, :)
   n_mol = length(sys) ÷ 3
   for m in 1:n_mol
      O  = X[3m - 2]
      H1 = X[3m - 1]
      H2 = X[3m    ]
      rOH1 = norm(H1 - O)
      rOH2 = norm(H2 - O)
      @test ustrip(u"Å", rOH1) ≈ 0.9572 atol = 1e-10
      @test ustrip(u"Å", rOH2) ≈ 0.9572 atol = 1e-10
      cosθ = (H1 - O) ⋅ (H2 - O) / (rOH1 * rOH2)
      @test acos(cosθ) ≈ 104.52 * π / 180 atol = 1e-10
   end
end

@testset "tip3p_water — Poisson-disk respects d_min (minimum-image PBC)" begin
   box   = 15.0u"Å"
   d_min = 2.7u"Å"
   sys   = tip3p_water(box; d_min = d_min, rng = MersenneTwister(0))
   X     = position(sys, :)
   # Oxygens are at positions 1, 4, 7, ...
   O_idx = 1:3:length(sys)
   bv    = ustrip(u"Å", box)
   for (a_idx, i) in enumerate(O_idx), j in O_idx[a_idx+1:end]
      dx_unitful = X[i] - X[j]
      # minimum-image
      dx = [ustrip(u"Å", dx_unitful[k]) - bv * round(ustrip(u"Å", dx_unitful[k]) / bv) for k in 1:3]
      @test sqrt(sum(abs2, dx)) ≥ ustrip(u"Å", d_min) - 1e-9
   end
end

@testset "tip3p_water — rng reproducibility" begin
   a = tip3p_water(10.0u"Å"; rng = MersenneTwister(11))
   b = tip3p_water(10.0u"Å"; rng = MersenneTwister(11))
   @test position(a, :) == position(b, :)

   c = tip3p_water(10.0u"Å"; rng = MersenneTwister(12))
   @test position(a, :) != position(c, :)
end

@testset "tip3p_water — bare-Real arguments interpreted as Å" begin
   a = tip3p_water(10.0u"Å"; rng = MersenneTwister(5))
   b = tip3p_water(10.0;     rng = MersenneTwister(5))
   @test position(a, :) == position(b, :)
end

@testset "tip3p_water — custom charge_label" begin
   sys = tip3p_water(10.0u"Å"; charge_label = :q, rng = MersenneTwister(0))
   @test sys[1].data[:q] == -0.834u"e_au"
   @test sys[2].data[:q] == +0.417u"e_au"
end
